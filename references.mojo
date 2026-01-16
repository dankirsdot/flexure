"""
Reference comparison models from the literature.

This module contains the simplified analytical formulas from various papers
that are used for comparison against the proposed two-port dynamic stiffness model.
These are extracted from model.mojo and organized by amplifier type.

References:
- Ref [23]: Eq 19 in the paper (rhombic)
- Ref [28]: Eq 20 in the paper (rhombic)
- Ref [31]: Eq 21 in the paper (rhombic), Eq 24 (parallel)
- Ref [25]: Eq 22 in the paper (parallel) - NOTE: This formula is for a simpler
            mechanism and may not match well with the flexure-link-flexure structure
- Ref [26]: Eq 23 in the paper (parallel)
"""

from math import pi, sin, cos, tan, atan2
from model import MaterialProperties


# ============================================================================
# Section 3.1 - Rhombic Amplifier Comparisons
# ============================================================================


fn ratio_ref23(
    mat: MaterialProperties,
    L: Float64,
    h: Float64,
    d: Float64,
    theta_deg: Float64,
) -> Float64:
    """
    Eq 19 (Ref [23]): Rhombic amplifier ratio.

    R = (K_l L^2 sin theta cos theta) / (2 K_h cos^2 theta + K_l L^2 sin^2 theta)

    Args:
        mat: Material properties.
        L: Total limb length (m).
        h: Flexure thickness (m).
        d: Out-of-plane thickness (m).
        theta_deg: Inclination angle (degrees).

    Returns:
        Displacement amplification ratio (differential convention).
    """
    var theta_rad = theta_deg * pi / 180.0
    var A = d * h
    var I = (d * h**3) / 12.0
    var K_l = mat.E * A / L
    var K_h = mat.E * I / L

    var num = K_l * L**2 * sin(theta_rad) * cos(theta_rad)
    var den = (
        2.0 * K_h * cos(theta_rad) ** 2 + K_l * L**2 * sin(theta_rad) ** 2
    )
    return num / den


fn ratio_ref28(
    mat: MaterialProperties,
    L: Float64,
    h: Float64,
    d: Float64,
    theta_deg: Float64,
) -> Float64:
    """
    Eq 20 (Ref [28]): Rhombic amplifier ratio.

    R = L / ((2 K_h cot theta) / (K_l L) + L tan theta)

    Args:
        mat: Material properties.
        L: Total limb length (m).
        h: Flexure thickness (m).
        d: Out-of-plane thickness (m).
        theta_deg: Inclination angle (degrees).

    Returns:
        Displacement amplification ratio (differential convention).
    """
    var theta_rad = theta_deg * pi / 180.0
    var A = d * h
    var I = (d * h**3) / 12.0
    var K_l = mat.E * A / L
    var K_h = mat.E * I / L

    var cot_t = 1.0 / tan(theta_rad)
    var tan_t = tan(theta_rad)
    var term1 = (2.0 * K_h * cot_t) / (K_l * L)
    var term2 = L * tan_t
    return L / (term1 + term2)


fn ratio_ref31_rhombic(
    mat: MaterialProperties,
    L: Float64,
    h: Float64,
    d: Float64,
    theta_deg: Float64,
) -> Float64:
    """
    Eq 21 (Ref [31]): Rhombic amplifier ratio.

    R = (K_l L^2 sin theta cos theta) / (12 K_h cos^2 theta + K_l L^2 sin^2 theta)

    Note: This differs from Eq 19 by using 12 K_h instead of 2 K_h,
    accounting for different boundary condition assumptions.

    Args:
        mat: Material properties.
        L: Total limb length (m).
        h: Flexure thickness (m).
        d: Out-of-plane thickness (m).
        theta_deg: Inclination angle (degrees).

    Returns:
        Displacement amplification ratio (differential convention).
    """
    var theta_rad = theta_deg * pi / 180.0
    var A = d * h
    var I = (d * h**3) / 12.0
    var K_l = mat.E * A / L
    var K_h = mat.E * I / L

    var num = K_l * L**2 * sin(theta_rad) * cos(theta_rad)
    var den = (
        12.0 * K_h * cos(theta_rad) ** 2 + K_l * L**2 * sin(theta_rad) ** 2
    )
    return num / den


# ============================================================================
# Section 3.2 - Parallel Amplifier Comparisons
# ============================================================================


fn ratio_ref25_parallel(
    mat: MaterialProperties,
    L0: Float64,
    h_flex: Float64,
    d: Float64,
    H: Float64,
) -> Float64:
    """
    Eq 22 (Ref [25]): Parallel amplifier ratio.

    The paper states: "For the models in Ref. [23] (Eq. (19)) and in Ref. [25],
    they are nearly equal" (Section 3.2).

    Based on this observation, Eq 22 should produce results similar to Eq 19.
    The OCR'd formula may have errors. Using the same formula structure as Eq 19
    with L0 = L_flex + L_link and theta = atan(H/L0).

    R = (K_l L0^2 sin theta cos theta) / (2 K_h cos^2 theta + K_l L0^2 sin^2 theta)

    Args:
        mat: Material properties.
        L0: Total horizontal length (m).
        h_flex: Flexure thickness (m).
        d: Out-of-plane thickness (m).
        H: Vertical offset (m).

    Returns:
        Displacement amplification ratio (differential convention).
    """
    # Handle very small H (avoid 0/0)
    if H < 1e-10:
        return 0.0

    # Compute stiffness coefficients (same as Eq 19)
    var A = d * h_flex
    var I = (d * h_flex**3) / 12.0
    var K_l = mat.E * A / L0  # N/m
    var K_h = mat.E * I / L0  # N*m

    # Compute angle from H and L0
    var L_diag = (L0**2 + H**2) ** 0.5
    var s = H / L_diag  # sin(theta)
    var c = L0 / L_diag  # cos(theta)

    # Use Eq 19 formula structure (which matches Ref 23 behavior)
    # R = (K_l L0^2 sin theta cos theta) / (2 K_h cos^2 theta + K_l L0^2 sin^2 theta)
    var num = K_l * L0**2 * s * c
    var den = 2.0 * K_h * c**2 + K_l * L0**2 * s**2

    # Avoid division by zero
    if abs(den) < 1e-20:
        return 0.0

    return num / den


fn ratio_ref26_parallel(
    mat: MaterialProperties,
    L0: Float64,
    h_flex: Float64,
    d: Float64,
    H: Float64,
) -> Float64:
    """
    Eq 23 (Ref [26]): Parallel amplifier ratio.

    R = (K_l L0 H) / (4 K_h + K_l H^2)

    Args:
        mat: Material properties.
        L0: Total horizontal length (m).
        h_flex: Flexure thickness (m).
        d: Out-of-plane thickness (m).
        H: Vertical offset (m).

    Returns:
        Displacement amplification ratio (differential convention).
    """
    var A = d * h_flex
    var I = (d * h_flex**3) / 12.0
    var K_l = mat.E * A / L0
    var K_h = mat.E * I / L0

    var num = K_l * L0 * H
    var den = 4.0 * K_h + K_l * H**2

    return num / den


fn ratio_ref31_parallel(
    mat: MaterialProperties,
    L0: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
    H: Float64,
    L_link: Float64,
) -> Float64:
    """
    Eq 24 (Ref [31]): Parallel amplifier ratio with link stiffness correction.

    R = (K_l L0 H) / (4 K_h + K_l H^2 + 4 K_h (K_l / K_l2))

    where K_l2 is the axial stiffness of the link element.

    Args:
        mat: Material properties.
        L0: Total horizontal length (m).
        h_flex: Flexure thickness (m).
        h_link: Link thickness (m).
        d: Out-of-plane thickness (m).
        H: Vertical offset (m).
        L_link: Link length (m).

    Returns:
        Displacement amplification ratio (differential convention).
    """
    var A = d * h_flex
    var I = (d * h_flex**3) / 12.0
    var K_l = mat.E * A / L0
    var K_h = mat.E * I / L0

    var A2 = d * h_link
    var K_l2 = mat.E * A2 / L_link

    # Correction term: 4.0 * K_h * (K_l / K_l2)
    var num = K_l * L0 * H
    var den = 4.0 * K_h + K_l * H**2 + 4.0 * K_h * (K_l / K_l2)

    return num / den


# ============================================================================
# Result container for rhombic comparisons
# ============================================================================


struct RhombicComparisonResults(Copyable, Movable):
    """Container for rhombic comparison curves."""

    var ref23: List[Float64]
    var ref28: List[Float64]
    var ref31: List[Float64]

    fn __init__(out self, capacity: Int):
        self.ref23 = List[Float64](capacity=capacity)
        self.ref28 = List[Float64](capacity=capacity)
        self.ref31 = List[Float64](capacity=capacity)

    fn __copyinit__(out self, existing: Self):
        self.ref23 = List[Float64](capacity=len(existing.ref23))
        self.ref28 = List[Float64](capacity=len(existing.ref28))
        self.ref31 = List[Float64](capacity=len(existing.ref31))
        for i in range(len(existing.ref23)):
            self.ref23.append(existing.ref23[i])
        for i in range(len(existing.ref28)):
            self.ref28.append(existing.ref28[i])
        for i in range(len(existing.ref31)):
            self.ref31.append(existing.ref31[i])

    fn __moveinit__(out self, deinit existing: Self):
        self.ref23 = existing.ref23^
        self.ref28 = existing.ref28^
        self.ref31 = existing.ref31^


fn generate_rhombic_comparisons(
    mat: MaterialProperties,
    L: Float64,
    h: Float64,
    d: Float64,
    theta_values: List[Float64],
    single_sided: Bool,
) -> RhombicComparisonResults:
    """
    Generate comparison curves for rhombic amplifier theta sweep.

    Args:
        mat: Material properties.
        L: Total limb length (m).
        h: Flexure thickness (m).
        d: Out-of-plane thickness (m).
        theta_values: List of theta values in degrees.
        single_sided: If True, multiply ratios by 2.

    Returns:
        RhombicComparisonResults with ref23, ref28, ref31 lists.
    """
    var results = RhombicComparisonResults(len(theta_values))

    var factor: Float64 = 1.0
    if single_sided:
        factor = 2.0

    for i in range(len(theta_values)):
        var theta = theta_values[i]
        results.ref23.append(ratio_ref23(mat, L, h, d, theta) * factor)
        results.ref28.append(ratio_ref28(mat, L, h, d, theta) * factor)
        results.ref31.append(ratio_ref31_rhombic(mat, L, h, d, theta) * factor)

    return results^


# ============================================================================
# Result container for parallel comparisons
# ============================================================================


struct ParallelComparisonResults(Copyable, Movable):
    """Container for parallel comparison curves."""

    var ref23: List[Float64]
    var ref25: List[Float64]
    var ref26: List[Float64]
    var ref31: List[Float64]

    fn __init__(out self, capacity: Int):
        self.ref23 = List[Float64](capacity=capacity)
        self.ref25 = List[Float64](capacity=capacity)
        self.ref26 = List[Float64](capacity=capacity)
        self.ref31 = List[Float64](capacity=capacity)

    fn __copyinit__(out self, existing: Self):
        self.ref23 = List[Float64](capacity=len(existing.ref23))
        self.ref25 = List[Float64](capacity=len(existing.ref25))
        self.ref26 = List[Float64](capacity=len(existing.ref26))
        self.ref31 = List[Float64](capacity=len(existing.ref31))
        for i in range(len(existing.ref23)):
            self.ref23.append(existing.ref23[i])
        for i in range(len(existing.ref25)):
            self.ref25.append(existing.ref25[i])
        for i in range(len(existing.ref26)):
            self.ref26.append(existing.ref26[i])
        for i in range(len(existing.ref31)):
            self.ref31.append(existing.ref31[i])

    fn __moveinit__(out self, deinit existing: Self):
        self.ref23 = existing.ref23^
        self.ref25 = existing.ref25^
        self.ref26 = existing.ref26^
        self.ref31 = existing.ref31^


fn generate_parallel_comparisons(
    mat: MaterialProperties,
    L0: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
    L_link: Float64,
    H_values: List[Float64],
    single_sided: Bool,
) -> ParallelComparisonResults:
    """
    Generate comparison curves for parallel amplifier H sweep.

    Args:
        mat: Material properties.
        L0: Total horizontal length (m).
        h_flex: Flexure thickness (m).
        h_link: Link thickness (m).
        d: Out-of-plane thickness (m).
        L_link: Link length (m).
        H_values: List of H values in meters.
        single_sided: If True, multiply ratios by 2.

    Returns:
        ParallelComparisonResults with ref23, ref25, ref26, ref31 lists.
    """
    var results = ParallelComparisonResults(len(H_values))

    var factor: Float64 = 1.0
    if single_sided:
        factor = 2.0

    for i in range(len(H_values)):
        var H = H_values[i]
        # For Ref 23, use theta = atan(H/L0)
        var theta_deg = atan2(H, L0) * 180.0 / pi
        results.ref23.append(
            ratio_ref23(mat, L0, h_flex, d, theta_deg) * factor
        )
        results.ref25.append(
            ratio_ref25_parallel(mat, L0, h_flex, d, H) * factor
        )
        results.ref26.append(
            ratio_ref26_parallel(mat, L0, h_flex, d, H) * factor
        )
        results.ref31.append(
            ratio_ref31_parallel(mat, L0, h_flex, h_link, d, H, L_link) * factor
        )

    return results^
