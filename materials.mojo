"""
Material properties library for flexure amplifier analysis.

This module provides commonly used materials as pre-defined constants,
avoiding magic numbers in user code and ensuring consistency across analyses.

Usage:
    from materials import aluminum_6061, steel_304ss
    
    let mat = aluminum_6061()
    let results = compute_results(mat, geom, mass)
"""

from model import MaterialProperties


# ============================================================================
# Common Engineering Materials (as factory functions)
# ============================================================================


fn aluminum_6061() -> MaterialProperties:
    """6061-T6 Aluminum: E=71 GPa, ρ=2770 kg/m^3. Default material in the paper.
    """
    return MaterialProperties(E=71.0e9, rho=2770.0)


fn aluminum_7075() -> MaterialProperties:
    """7075-T6 Aluminum: E=71.7 GPa, ρ=2810 kg/m^3. Higher strength variant."""
    return MaterialProperties(E=71.7e9, rho=2810.0)


fn steel_304ss() -> MaterialProperties:
    """304 Stainless Steel: E=193 GPa, ρ=8000 kg/m^3."""
    return MaterialProperties(E=193.0e9, rho=8000.0)


fn steel_17_4ph() -> MaterialProperties:
    """17-4 PH Stainless Steel: E=197 GPa, ρ=7800 kg/m^3. Common for flexures.
    """
    return MaterialProperties(E=197.0e9, rho=7800.0)


fn titanium_ti6al4v() -> MaterialProperties:
    """Ti-6Al-4V Titanium: E=113.8 GPa, ρ=4430 kg/m^3. High strength-to-weight.
    """
    return MaterialProperties(E=113.8e9, rho=4430.0)


fn beryllium_copper() -> MaterialProperties:
    """Beryllium Copper (C17200): E=131 GPa, ρ=8250 kg/m^3. Fatigue resistant.
    """
    return MaterialProperties(E=131.0e9, rho=8250.0)


fn invar_36() -> MaterialProperties:
    """Invar 36: E=141 GPa, ρ=8050 kg/m^3. Low thermal expansion."""
    return MaterialProperties(E=141.0e9, rho=8050.0)


# ============================================================================
# Factory Function for Custom Materials
# ============================================================================


fn custom_material(E: Float64, rho: Float64) -> MaterialProperties:
    """
    Create a custom material with given properties.

    Args:
        E: Young's modulus in Pascals (Pa).
        rho: Density in kg/m^3.

    Returns:
        MaterialProperties instance ready for analysis.

    Example:
        let polymer = custom_material(E=3.0e9, rho=1200.0)
    """
    return MaterialProperties(E=E, rho=rho)
