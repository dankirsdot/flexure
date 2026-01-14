"""
Generic parameter sweep utilities for flexure amplifiers.

This module provides reusable sweep functions that eliminate the duplicated
sweep code in the original main_modular.mojo. Uses configuration structs
that capture the sweep parameters.

Usage:
    from sweeps import (
        RhombicSweepConfig, ParallelSweepConfig, AlignedSweepConfig,
        sweep_rhombic_theta, sweep_parallel_H, sweep_aligned_theta,
        Sweep1DResults
    )
    
    let config = make_rhombic_config(L=15e-3, h=1e-3, d=10e-3, L_flex=3e-3)
    let results = sweep_rhombic_theta(aluminum_6061(), config, mass, 0.01, 30.0, 500)
"""

from model import (
    MaterialProperties,
    GeometricParameters,
    calculate_amplification_ratio,
    calculate_input_stiffness,
    calculate_natural_frequencies,
)
from geometry import (
    rhombic, parallel, aligned,
    InputPortParams, default_input_port,
)
from analysis import OutputMass


# ============================================================================
# Sweep Results Container
# ============================================================================

struct Sweep1DResults(Copyable, Movable):
    """
    Results from a 1D parameter sweep.
    
    Contains all computed metrics at each parameter value, ready for plotting.
    
    Attributes:
        param_values: The swept parameter values.
        param_name: Name for axis labels (e.g., "theta", "H").
        param_unit: Unit for axis labels (e.g., "deg", "mm").
        ratios: Amplification ratio at each parameter value.
        stiffnesses: Input stiffness in N/um (already scaled from N/m).
        f1, f2, f3: Natural frequencies at each parameter value.
    """
    var param_values: List[Float64]
    var param_name: String
    var param_unit: String
    var ratios: List[Float64]
    var stiffnesses: List[Float64]
    var f1: List[Float64]
    var f2: List[Float64]
    var f3: List[Float64]
    
    fn __init__(out self, capacity: Int, param_name: String, param_unit: String):
        """
        Create sweep results with pre-allocated capacity.
        
        Args:
            capacity: Expected number of sweep points.
            param_name: Name of the swept parameter.
            param_unit: Unit of the swept parameter.
        """
        self.param_values = List[Float64](capacity=capacity)
        self.param_name = param_name
        self.param_unit = param_unit
        self.ratios = List[Float64](capacity=capacity)
        self.stiffnesses = List[Float64](capacity=capacity)
        self.f1 = List[Float64](capacity=capacity)
        self.f2 = List[Float64](capacity=capacity)
        self.f3 = List[Float64](capacity=capacity)
    
    fn __copyinit__(out self, existing: Self):
        """Copy constructor."""
        self.param_values = List[Float64](capacity=len(existing.param_values))
        self.param_name = existing.param_name
        self.param_unit = existing.param_unit
        self.ratios = List[Float64](capacity=len(existing.ratios))
        self.stiffnesses = List[Float64](capacity=len(existing.stiffnesses))
        self.f1 = List[Float64](capacity=len(existing.f1))
        self.f2 = List[Float64](capacity=len(existing.f2))
        self.f3 = List[Float64](capacity=len(existing.f3))
        
        for i in range(len(existing.param_values)):
            self.param_values.append(existing.param_values[i])
        for i in range(len(existing.ratios)):
            self.ratios.append(existing.ratios[i])
        for i in range(len(existing.stiffnesses)):
            self.stiffnesses.append(existing.stiffnesses[i])
        for i in range(len(existing.f1)):
            self.f1.append(existing.f1[i])
        for i in range(len(existing.f2)):
            self.f2.append(existing.f2[i])
        for i in range(len(existing.f3)):
            self.f3.append(existing.f3[i])
    
    fn __moveinit__(out self, deinit existing: Self):
        """Move constructor."""
        self.param_values = existing.param_values^
        self.param_name = existing.param_name^
        self.param_unit = existing.param_unit^
        self.ratios = existing.ratios^
        self.stiffnesses = existing.stiffnesses^
        self.f1 = existing.f1^
        self.f2 = existing.f2^
        self.f3 = existing.f3^
    
    fn add_point(
        mut self,
        param_value: Float64,
        ratio: Float64,
        stiffness_N_per_m: Float64,
        var freqs: List[Float64],
    ):
        """
        Add a single sweep point to the results.
        
        Args:
            param_value: The parameter value at this point.
            ratio: Amplification ratio.
            stiffness_N_per_m: Input stiffness in N/m (will be converted to N/um).
            freqs: List of natural frequencies.
        """
        self.param_values.append(param_value)
        self.ratios.append(ratio)
        self.stiffnesses.append(stiffness_N_per_m / 1.0e6) # Convert to N/um
        
        if len(freqs) >= 1:
            self.f1.append(freqs[0])
        else:
            self.f1.append(0.0)
        
        if len(freqs) >= 2:
            self.f2.append(freqs[1])
        else:
            self.f2.append(0.0)
        
        if len(freqs) >= 3:
            self.f3.append(freqs[2])
        else:
            self.f3.append(0.0)


# ============================================================================
# Rhombic Amplifier Sweep
# ============================================================================

@fieldwise_init
struct RhombicSweepConfig(Copyable, Movable):
    """
    Configuration for rhombic amplifier geometry.
    
    Captures all parameters except the one being swept (theta).
    
    Attributes:
        L: Total limb length (m).
        h: Flexure thickness (m).
        d: Out-of-plane depth (m).
        L_flex: Flexure hinge length (m).
        input_port: Input port dimensions.
    """
    var L: Float64
    var h: Float64
    var d: Float64
    var L_flex: Float64
    var input_port: InputPortParams


fn make_rhombic_config(
    L: Float64,
    h: Float64,
    d: Float64,
    L_flex: Float64,
) -> RhombicSweepConfig:
    """Create rhombic sweep config with default input port."""
    return RhombicSweepConfig(
        L=L, h=h, d=d, L_flex=L_flex,
        input_port=default_input_port()
    )


fn sweep_rhombic_theta(
    mat: MaterialProperties,
    config: RhombicSweepConfig,
    output_mass: OutputMass,
    theta_start: Float64,
    theta_end: Float64,
    num_steps: Int,
    single_sided: Bool = True,
    input_mass: Float64 = 0.0,
    verbose: Bool = True,
) -> Sweep1DResults:
    """
    Sweep theta (inclination angle) for a rhombic amplifier.
    
    Reproduces Figure 4 from the paper when used with appropriate parameters.
    
    Args:
        mat: Material properties.
        config: Rhombic geometry configuration (all params except theta).
        output_mass: Output port mass and inertia.
        theta_start: Starting angle (degrees).
        theta_end: Ending angle (degrees).
        num_steps: Number of sweep points.
        single_sided: If True, use single-sided convention (x2).
        input_mass: Input port mass (kg). Default: 0.
        verbose: If True, print progress every 20 steps.
    
    Returns:
        Sweep1DResults with all computed metrics.
    """
    var results = Sweep1DResults(num_steps, "theta", "deg")
    var step_size = (theta_end - theta_start) / Float64(num_steps - 1)
    
    for i in range(num_steps):
        var theta_deg = theta_start + Float64(i) * step_size
        
        if verbose and i % 20 == 0:
            print("Rhombic sweep: step", i, "/", num_steps)
        
        # Generate geometry at this theta
        var geom = rhombic(
            L=config.L,
            h=config.h,
            d=config.d,
            theta_deg=theta_deg,
            L_flex=config.L_flex,
            input_port=config.input_port,
        )
        
        # Compute metrics
        var ratio = calculate_amplification_ratio(mat, geom, single_sided)
        var stiffness = calculate_input_stiffness(mat, geom, single_sided)
        var freqs = calculate_natural_frequencies(
            mat, geom, output_mass.m, output_mass.J, input_mass
        )
        
        results.add_point(theta_deg, ratio, stiffness, freqs^)
    
    return results^


# ============================================================================
# Parallel Amplifier Sweep
# ============================================================================

@fieldwise_init
struct ParallelSweepConfig(Copyable, Movable):
    """
    Configuration for parallel amplifier geometry.
    
    Captures all parameters except the one being swept (H).
    
    Attributes:
        L_flex: Horizontal flexure length (m).
        L_link: Link horizontal projection (m).
        h_flex: Flexure thickness (m).
        h_link: Link thickness (m).
        d: Out-of-plane depth (m).
        input_port: Input port dimensions.
    """
    var L_flex: Float64
    var L_link: Float64
    var h_flex: Float64
    var h_link: Float64
    var d: Float64
    var input_port: InputPortParams


fn make_parallel_config(
    L_flex: Float64,
    L_link: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
) -> ParallelSweepConfig:
    """Create parallel sweep config with default input port."""
    return ParallelSweepConfig(
        L_flex=L_flex, L_link=L_link, h_flex=h_flex, h_link=h_link, d=d,
        input_port=default_input_port()
    )


fn sweep_parallel_H(
    mat: MaterialProperties,
    config: ParallelSweepConfig,
    output_mass: OutputMass,
    H_start: Float64,
    H_end: Float64,
    num_steps: Int,
    single_sided: Bool = True,
    input_mass: Float64 = 0.0,
    verbose: Bool = True,
) -> Sweep1DResults:
    """
    Sweep H (vertical offset) for a parallel amplifier.
    
    Reproduces Figure 5 from the paper when used with appropriate parameters.
    Note: param_values are stored in mm for plotting convenience.
    
    Args:
        mat: Material properties.
        config: Parallel geometry configuration (all params except H).
        output_mass: Output port mass and inertia.
        H_start: Starting offset (m).
        H_end: Ending offset (m).
        num_steps: Number of sweep points.
        single_sided: If True, use single-sided convention (x2).
        input_mass: Input port mass (kg). Default: 0.
        verbose: If True, print progress every 20 steps.
    
    Returns:
        Sweep1DResults with param_values in mm.
    """
    var results = Sweep1DResults(num_steps, "H", "mm")
    var step_size = (H_end - H_start) / Float64(num_steps - 1)
    
    for i in range(num_steps):
        var H = H_start + Float64(i) * step_size
        
        if verbose and i % 20 == 0:
            print("Parallel sweep: step", i, "/", num_steps)
        
        # Generate geometry at this H
        var geom = parallel(
            L_flex=config.L_flex,
            L_link=config.L_link,
            h_flex=config.h_flex,
            h_link=config.h_link,
            d=config.d,
            H=H,
            input_port=config.input_port,
        )
        
        # Compute metrics
        var ratio = calculate_amplification_ratio(mat, geom, single_sided)
        var stiffness = calculate_input_stiffness(mat, geom, single_sided)
        var freqs = calculate_natural_frequencies(
            mat, geom, output_mass.m, output_mass.J, input_mass
        )
        
        # Store H in mm for plotting
        results.add_point(H * 1000.0, ratio, stiffness, freqs^)
    
    return results^


# ============================================================================
# Aligned Amplifier Sweep
# ============================================================================

@fieldwise_init
struct AlignedSweepConfig(Copyable, Movable):
    """
    Configuration for aligned amplifier geometry.
    
    Captures all parameters except the one being swept (theta).
    
    Attributes:
        L_flex: Flexure length (m).
        L_link: Central link length (m).
        h_flex: Flexure thickness (m).
        h_link: Link thickness (m).
        d: Out-of-plane depth (m).
        input_port: Input port dimensions.
    """
    var L_flex: Float64
    var L_link: Float64
    var h_flex: Float64
    var h_link: Float64
    var d: Float64
    var input_port: InputPortParams


fn make_aligned_config(
    L_flex: Float64,
    L_link: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
) -> AlignedSweepConfig:
    """Create aligned sweep config with default input port."""
    return AlignedSweepConfig(
        L_flex=L_flex, L_link=L_link, h_flex=h_flex, h_link=h_link, d=d,
        input_port=default_input_port()
    )


fn sweep_aligned_theta(
    mat: MaterialProperties,
    config: AlignedSweepConfig,
    output_mass: OutputMass,
    theta_start: Float64,
    theta_end: Float64,
    num_steps: Int,
    single_sided: Bool = True,
    input_mass: Float64 = 0.0,
    verbose: Bool = True,
) -> Sweep1DResults:
    """
    Sweep theta (inclination angle) for an aligned amplifier.
    
    Reproduces Figure 6 from the paper when used with appropriate parameters.
    
    Args:
        mat: Material properties.
        config: Aligned geometry configuration (all params except theta).
        output_mass: Output port mass and inertia.
        theta_start: Starting angle (degrees).
        theta_end: Ending angle (degrees).
        num_steps: Number of sweep points.
        single_sided: If True, use single-sided convention (x2).
        input_mass: Input port mass (kg). Default: 0.
        verbose: If True, print progress every 20 steps.
    
    Returns:
        Sweep1DResults with all computed metrics.
    """
    var results = Sweep1DResults(num_steps, "theta", "deg")
    var step_size = (theta_end - theta_start) / Float64(num_steps - 1)
    
    for i in range(num_steps):
        var theta_deg = theta_start + Float64(i) * step_size
        
        if verbose and i % 20 == 0:
            print("Aligned sweep: step", i, "/", num_steps)
        
        # Generate geometry at this theta
        var geom = aligned(
            L_flex=config.L_flex,
            L_link=config.L_link,
            h_flex=config.h_flex,
            h_link=config.h_link,
            d=config.d,
            theta_deg=theta_deg,
            input_port=config.input_port,
        )
        
        # Compute metrics
        var ratio = calculate_amplification_ratio(mat, geom, single_sided)
        var stiffness = calculate_input_stiffness(mat, geom, single_sided)
        var freqs = calculate_natural_frequencies(
            mat, geom, output_mass.m, output_mass.J, input_mass
        )
        
        results.add_point(theta_deg, ratio, stiffness, freqs^)
    
    return results^


# ============================================================================
# Utility: Extract values from sweep results (for comparison curves)
# ============================================================================

fn get_theta_values(results: Sweep1DResults) -> List[Float64]:
    """
    Extract parameter values as a new list (for use with reference functions).
    
    Args:
        results: Sweep results from a theta sweep.
    
    Returns:
        Copy of param_values list.
    """
    var values = List[Float64](capacity=len(results.param_values))
    for i in range(len(results.param_values)):
        values.append(results.param_values[i])
    return values^


fn get_H_values_meters(results: Sweep1DResults) -> List[Float64]:
    """
    Extract H values in meters (sweep stores in mm).
    
    Args:
        results: Sweep results from an H sweep.
    
    Returns:
        List of H values in meters.
    """
    var values = List[Float64](capacity=len(results.param_values))
    for i in range(len(results.param_values)):
        values.append(results.param_values[i] / 1000.0) # mm -> m
    return values^
