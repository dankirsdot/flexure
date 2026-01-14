"""
Unified geometry construction for bridge-type flexure amplifiers.

This module embodies the paper's core insight: all three amplifier types
(rhombic, parallel, aligned) share the same 4-beam topology and can be
described by a single GeometricParameters struct. The factory functions
here translate high-level design parameters into the universal representation.

Usage:
    from geometry import rhombic, parallel, aligned, InputPortParams
    
    # Rhombic amplifier at 10 degrees
    let geom = rhombic(L=15e-3, h=1e-3, d=10e-3, theta_deg=10.0, L_flex=3e-3)
    
    # Parallel amplifier with 2mm offset
    let geom = parallel(L_flex=3e-3, L_link=9e-3, h_flex=0.5e-3, h_link=4e-3, d=10e-3, H=2e-3)
"""

from math import pi, sqrt, atan2
from model import GeometricParameters


# ============================================================================
# Input Port Configuration
# ============================================================================

@fieldwise_init
struct InputPortParams(Copyable, Movable):
    """
    Parameters for the input port rigid block.
    
    The input port is modeled as a short, stiff beam connecting the
    actuation point to the first flexure. Default values match the paper.
    
    Attributes:
        L0: Length of input block in x-direction (m).
        h0: Height of input block in y-direction (m).
    """
    var L0: Float64
    var h0: Float64


fn default_input_port() -> InputPortParams:
    """Default input port dimensions matching the paper (10mm x 5mm)."""
    return InputPortParams(L0=10.0e-3, h0=5.0e-3)


# ============================================================================
# Internal Helper
# ============================================================================

struct BeamParams:
    """Helper struct to return beam parameters."""
    var L: Float64
    var h: Float64
    var theta: Float64
    
    fn __init__(out self, L: Float64, h: Float64, theta: Float64):
        self.L = L
        self.h = h
        self.theta = theta


fn _compute_input_beam(input_port: InputPortParams) -> BeamParams:
    """
    Compute Beam 1 parameters from input port dimensions.
    
    The input port half is modeled as a diagonal beam from the actuation
    point to the corner where the limb begins.
    
    Returns:
        BeamParams with L_1, h_1, theta_1 (theta in radians).
    """
    var L_1 = sqrt((input_port.L0 / 2.0)**2 + (input_port.h0 / 2.0)**2)
    var h_1 = input_port.h0 # Treat as stiff short beam
    var theta_1 = atan2(input_port.L0 / 2.0, input_port.h0 / 2.0)
    return BeamParams(L_1, h_1, theta_1)


# ============================================================================
# Rhombic Amplifier (Section 3.1)
# ============================================================================

fn rhombic(
    L: Float64,
    h: Float64,
    d: Float64,
    theta_deg: Float64,
    L_flex: Float64,
    input_port: InputPortParams,
) -> GeometricParameters:
    """
    Rhombic-type amplifier geometry (Paper Section 3.1).
    
    The limb structure is:
        Input_Block -> Flexure1 -> Rigid_Link -> Flexure2 -> Output
    
    All three active beams (Flexure1, Rigid_Link, Flexure2) share the same
    inclination angle theta and thickness h. The rigid link length is computed
    as L_rigid = L - 2xL_flex.
    
    This is the classic diamond-shaped displacement amplifier where the
    amplification ratio approaches cot(theta) for small angles.
    
    Args:
        L: Total limb length from input corner to output (m).
        h: Uniform thickness of flexure hinges and link (m).
        d: Out-of-plane depth/thickness (m).
        theta_deg: Inclination angle in degrees.
        L_flex: Length of each flexure hinge (m).
        input_port: Input port block dimensions.
    
    Returns:
        GeometricParameters ready for analysis.
    """
    var theta_rad = theta_deg * pi / 180.0
    var L_rigid = L - 2.0 * L_flex
    
    # Beam 1: Input port half
    var input_beam = _compute_input_beam(input_port)
    var L_1 = input_beam.L
    var h_1 = input_beam.h
    var theta_1 = input_beam.theta

    # Beam 2: Flexure 1
    var L_2 = L_flex
    var h_2 = h
    var theta_2 = theta_rad

    # Beam 3: Rigid link
    var L_3 = L_rigid
    var h_3 = h
    var theta_3 = theta_rad

    # Beam 4: Flexure 2
    var L_4 = L_flex
    var h_4 = h
    var theta_4 = theta_rad

    return GeometricParameters(
        L_1=L_1, L_2=L_2, L_3=L_3, L_4=L_4,
        h_1=h_1, h_2=h_2, h_3=h_3, h_4=h_4,
        theta_1=theta_1, theta_2=theta_2, theta_3=theta_3, theta_4=theta_4,
        H=0.0, d=d
    )


fn rhombic_default(
    L: Float64,
    h: Float64,
    d: Float64,
    theta_deg: Float64,
    L_flex: Float64,
) -> GeometricParameters:
    """Rhombic geometry with default input port (10mm x 5mm)."""
    return rhombic(L, h, d, theta_deg, L_flex, default_input_port())


# ============================================================================
# Parallel Amplifier (Section 3.2)
# ============================================================================

fn parallel(
    L_flex: Float64,
    L_link: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
    H: Float64,
    input_port: InputPortParams,
) -> GeometricParameters:
    """
    Parallel-type amplifier geometry (Paper Section 3.2).
    
    The limb structure is:
        Input_Block -> Flexure1(horiz) -> Link(angled) -> Flexure2(horiz) -> Output
    
    The horizontal flexures (theta=0) are connected by an angled link at
    theta_link = atan(H / L_link). This geometry provides good stiffness
    characteristics with the offset H controlling the amplification.
    
    Args:
        L_flex: Length of horizontal flexure hinges (m).
        L_link: Horizontal projection of the angled link (m).
        h_flex: Thickness of flexure hinges (m).
        h_link: Thickness of the angled link (m).
        d: Out-of-plane depth/thickness (m).
        H: Vertical offset between input and output (m).
        input_port: Input port block dimensions.
    
    Returns:
        GeometricParameters ready for analysis.
    """
    # Beam 1: Input port half
    var input_beam = _compute_input_beam(input_port)
    var L_1 = input_beam.L
    var h_1 = input_beam.h
    var theta_1 = input_beam.theta

    # Beam 2: Flexure 1 (horizontal)
    var L_2 = L_flex
    var h_2 = h_flex
    var theta_2 = 0.0

    # Beam 3: Angled link
    var L_3 = sqrt(L_link**2 + H**2)
    var h_3 = h_link
    var theta_3 = atan2(H, L_link)

    # Beam 4: Flexure 2 (horizontal)
    var L_4 = L_flex
    var h_4 = h_flex
    var theta_4 = 0.0

    return GeometricParameters(
        L_1=L_1, L_2=L_2, L_3=L_3, L_4=L_4,
        h_1=h_1, h_2=h_2, h_3=h_3, h_4=h_4,
        theta_1=theta_1, theta_2=theta_2, theta_3=theta_3, theta_4=theta_4,
        H=H, d=d
    )


fn parallel_default(
    L_flex: Float64,
    L_link: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
    H: Float64,
) -> GeometricParameters:
    """Parallel geometry with default input port (10mm x 5mm)."""
    return parallel(L_flex, L_link, h_flex, h_link, d, H, default_input_port())


# ============================================================================
# Aligned Amplifier (Section 3.3)
# ============================================================================

fn aligned(
    L_flex: Float64,
    L_link: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
    theta_deg: Float64,
    input_port: InputPortParams,
) -> GeometricParameters:
    """
    Aligned-type amplifier geometry (Paper Section 3.3).
    
    The limb structure is:
        Input_Block -> Flexure1 -> Link -> Flexure2 -> Output
    
    Similar to the rhombic type, all three active beams share the same
    inclination angle theta. However, the flexures and link have distinct
    dimensions (h_flex vs h_link, L_flex vs L_link), allowing more
    design freedom.
    
    Args:
        L_flex: Length of flexure hinges (m).
        L_link: Length of the central link (m).
        h_flex: Thickness of flexure hinges (m).
        h_link: Thickness of the central link (m).
        d: Out-of-plane depth/thickness (m).
        theta_deg: Inclination angle in degrees.
        input_port: Input port block dimensions.
    
    Returns:
        GeometricParameters ready for analysis.
    """
    var theta_rad = theta_deg * pi / 180.0
    
    # Beam 1: Input port half
    var input_beam = _compute_input_beam(input_port)
    var L_1 = input_beam.L
    var h_1 = input_beam.h
    var theta_1 = input_beam.theta

    # Beam 2: Flexure 1
    var L_2 = L_flex
    var h_2 = h_flex
    var theta_2 = theta_rad

    # Beam 3: Central link
    var L_3 = L_link
    var h_3 = h_link
    var theta_3 = theta_rad

    # Beam 4: Flexure 2
    var L_4 = L_flex
    var h_4 = h_flex
    var theta_4 = theta_rad

    return GeometricParameters(
        L_1=L_1, L_2=L_2, L_3=L_3, L_4=L_4,
        h_1=h_1, h_2=h_2, h_3=h_3, h_4=h_4,
        theta_1=theta_1, theta_2=theta_2, theta_3=theta_3, theta_4=theta_4,
        H=0.0, d=d
    )


fn aligned_default(
    L_flex: Float64,
    L_link: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
    theta_deg: Float64,
) -> GeometricParameters:
    """Aligned geometry with default input port (10mm x 5mm)."""
    return aligned(L_flex, L_link, h_flex, h_link, d, theta_deg, default_input_port())


# ============================================================================
# Custom Bridge Amplifier
# ============================================================================

fn bridge_amplifier(
    L_1: Float64, h_1: Float64, theta_1_deg: Float64,
    L_2: Float64, h_2: Float64, theta_2_deg: Float64,
    L_3: Float64, h_3: Float64, theta_3_deg: Float64,
    L_4: Float64, h_4: Float64, theta_4_deg: Float64,
    d: Float64,
) -> GeometricParameters:
    """
    Fully custom 4-beam bridge amplifier geometry.
    
    For advanced users who want complete control over each beam segment.
    This bypasses the high-level factories and directly specifies all
    12 beam parameters (4 beams x 3 parameters each).
    
    All angles are specified in degrees for user convenience.
    """
    var deg_to_rad = pi / 180.0
    
    return GeometricParameters(
        L_1=L_1, L_2=L_2, L_3=L_3, L_4=L_4,
        h_1=h_1, h_2=h_2, h_3=h_3, h_4=h_4,
        theta_1=theta_1_deg * deg_to_rad,
        theta_2=theta_2_deg * deg_to_rad,
        theta_3=theta_3_deg * deg_to_rad,
        theta_4=theta_4_deg * deg_to_rad,
        H=0.0, d=d
    )
