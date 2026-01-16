"""
Unified analysis functions for flexure amplifiers.

This module provides a single compute_results() function that replaces
the three separate solve_*_amplifier() functions in the original code.
It computes static (amplification ratio, input stiffness) and dynamic
(natural frequencies) properties for any amplifier geometry.

Usage:
    from analysis import compute_results, OutputMass, AmplifierResults
    from geometry import rhombic_default
    from materials import aluminum_6061
    
    let geom = rhombic_default(L=15e-3, h=1e-3, d=10e-3, theta_deg=10.0, L_flex=3e-3)
    let mass = make_output_mass_from_block(4e-3, 6e-3, 10e-3, aluminum_6061().rho)
    let results = compute_results(aluminum_6061(), geom, mass)
    results.dump()
"""

from model import (
    MaterialProperties,
    GeometricParameters,
    calculate_amplification_ratio,
    calculate_input_stiffness,
    calculate_natural_frequencies,
)


# ============================================================================
# Output Mass Configuration
# ============================================================================


@fieldwise_init
struct OutputMass(Copyable, Movable):
    """
    Lumped mass and rotational inertia at the output port.

    The output port typically has a rigid block (the "output shuttle")
    whose mass and inertia affect the natural frequencies but not the
    static amplification ratio.

    Attributes:
        m: Translational mass (kg).
        J: Rotational inertia about z-axis (kg*m^2).
    """

    var m: Float64
    var J: Float64


fn make_output_mass(m: Float64 = 0.0, J: Float64 = 0.0) -> OutputMass:
    """
    Create output mass parameters.

    Args:
        m: Translational mass (kg). Default: 0.
        J: Rotational inertia (kg*m^2). Default: 0.
    """
    return OutputMass(m=m, J=J)


fn make_output_mass_from_block(
    width: Float64, height: Float64, depth: Float64, rho: Float64
) -> OutputMass:
    """
    Compute mass and inertia from a rectangular block.

    Assumes a uniform density rectangular block with dimensions
    width x height x depth. The rotational inertia is computed
    about the centroid using the standard formula for a cuboid.

    Args:
        width: Block width in x-direction (m).
        height: Block height in y-direction (m).
        depth: Block depth in z-direction (m).
        rho: Material density (kg/m^3).

    Returns:
        OutputMass with computed m and J values.

    Example:
        # 4mm x 6mm x 10mm aluminum block
        let mass = make_output_mass_from_block(4e-3, 6e-3, 10e-3, 2770.0)
    """
    var volume = width * height * depth
    var m = volume * rho
    # Rotational inertia of cuboid about centroid: J = m/12 * (w^2 + h^2)
    var J = (m / 12.0) * (width**2 + height**2)
    return OutputMass(m=m, J=J)


# ============================================================================
# Analysis Results
# ============================================================================


struct AmplifierResults(Copyable, Movable):
    """
    Results from a single-point amplifier analysis.

    Contains the key performance metrics:
    - Amplification ratio (dimensionless)
    - Input stiffness (N/m)
    - First three natural frequencies (Hz)

    Attributes:
        amplification_ratio: Output displacement / input displacement.
        input_stiffness: Force required per unit input displacement (N/m).
        natural_frequencies: List of natural frequencies in Hz.
    """

    var amplification_ratio: Float64
    var input_stiffness: Float64
    var natural_frequencies: List[Float64]

    fn __init__(out self):
        """Create empty results (for initialization)."""
        self.amplification_ratio = 0.0
        self.input_stiffness = 0.0
        self.natural_frequencies = List[Float64]()

    fn __init__(
        out self, ratio: Float64, stiffness: Float64, var freqs: List[Float64]
    ):
        """
        Create results with specified values.

        Args:
            ratio: Amplification ratio.
            stiffness: Input stiffness (N/m).
            freqs: Natural frequencies (Hz).
        """
        self.amplification_ratio = ratio
        self.input_stiffness = stiffness
        self.natural_frequencies = freqs^

    fn __copyinit__(out self, existing: Self):
        """Copy constructor."""
        self.amplification_ratio = existing.amplification_ratio
        self.input_stiffness = existing.input_stiffness
        self.natural_frequencies = List[Float64](
            capacity=len(existing.natural_frequencies)
        )
        for i in range(len(existing.natural_frequencies)):
            self.natural_frequencies.append(existing.natural_frequencies[i])

    fn __moveinit__(out self, deinit existing: Self):
        """Move constructor."""
        self.amplification_ratio = existing.amplification_ratio
        self.input_stiffness = existing.input_stiffness
        self.natural_frequencies = existing.natural_frequencies^

    fn dump(self):
        """Print results to stdout in human-readable format."""
        print("=" * 50)
        print("Amplifier Analysis Results")
        print("=" * 50)
        print("Amplification Ratio:", self.amplification_ratio)
        print("Input Stiffness:", self.input_stiffness, "N/m")
        print("                ", self.input_stiffness / 1.0e6, "N/um")
        print("Natural Frequencies:")
        for i in range(len(self.natural_frequencies)):
            print("  Mode", i + 1, ":", self.natural_frequencies[i], "Hz")
        print("=" * 50)

    fn f1(self) -> Float64:
        """First natural frequency (Hz), or 0.0 if not computed."""
        if len(self.natural_frequencies) >= 1:
            return self.natural_frequencies[0]
        return 0.0

    fn f2(self) -> Float64:
        """Second natural frequency (Hz), or 0.0 if not computed."""
        if len(self.natural_frequencies) >= 2:
            return self.natural_frequencies[1]
        return 0.0

    fn f3(self) -> Float64:
        """Third natural frequency (Hz), or 0.0 if not computed."""
        if len(self.natural_frequencies) >= 3:
            return self.natural_frequencies[2]
        return 0.0


# ============================================================================
# Main Analysis Function
# ============================================================================


fn compute_results(
    mat: MaterialProperties,
    geom: GeometricParameters,
    output_mass: OutputMass,
    input_mass: Float64 = 0.0,
    single_sided: Bool = True,
) -> AmplifierResults:
    """
    Compute full static and dynamic analysis for an amplifier configuration.

    This is the unified analysis function that replaces the separate
    solve_rhombic_amplifier(), solve_parallel_amplifier(), and
    solve_aligned_amplifier() functions from the original code.

    Static Analysis:
    - Amplification ratio: Computed by applying equal and opposite forces
      at the input ports and measuring the output displacement.
    - Input stiffness: Force per unit input displacement.

    Dynamic Analysis:
    - Natural frequencies: Found by sweeping frequency and detecting
      zero crossings of the dynamic stiffness matrix determinant.

    Args:
        mat: Material properties (Young's modulus E, density rho).
        geom: Geometric parameters from geometry.mojo factory functions.
        output_mass: Lumped mass and inertia at the output port.
        input_mass: Lumped mass at each input port (kg). Default: 0.
        single_sided: If True, report single-sided ratio (multiply by 2).
            This is the default per the paper's convention where only
            one input port is actuated against a fixed reference.

    Returns:
        AmplifierResults containing ratio, stiffness, and frequencies.
    """
    # Static analysis
    var ratio = calculate_amplification_ratio(mat, geom, single_sided)
    var stiffness = calculate_input_stiffness(mat, geom, single_sided)

    # Dynamic analysis
    var freqs = calculate_natural_frequencies(
        mat, geom, output_mass.m, output_mass.J, input_mass
    )

    return AmplifierResults(ratio, stiffness, freqs^)
