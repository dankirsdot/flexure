"""
Example demonstrations for the flexure amplifier library.

This module contains example functions that reproduce the analyses from
the paper and demonstrate the library's capabilities. Each example is
designed to produce identical results to the original main_modular.mojo
for validation purposes.

Usage:
    from examples import (
        example_rhombic_single_point,
        example_rhombic_sweep,
        example_parallel_sweep,
        example_aligned_sweep,
    )
    
    # Run all examples
    example_rhombic_single_point()
    example_rhombic_sweep()
    example_parallel_sweep()
    example_aligned_sweep()
"""

from math import pi, tan

from materials import aluminum_6061
from geometry import rhombic_default, parallel_default, aligned_default
from analysis import compute_results, make_output_mass_from_block, OutputMass
from sweeps import (
    make_rhombic_config, make_parallel_config, make_aligned_config,
    sweep_rhombic_theta, sweep_parallel_H, sweep_aligned_theta,
    get_theta_values, get_H_values_meters,
)
from plotting import (
    plot_sweep_1d,
    plot_rhombic_with_comparisons,
    plot_parallel_with_comparisons,
)
from references import (
    generate_rhombic_comparisons,
    generate_parallel_comparisons,
)


# ============================================================================
# Common Parameters (matching main_modular.mojo)
# ============================================================================

fn _get_output_mass() -> OutputMass:
    """Standard output mass from the paper (4x6x10 mm aluminum block)."""
    var mat = aluminum_6061()
    return make_output_mass_from_block(4.0e-3, 6.0e-3, 10.0e-3, mat.rho)


# ============================================================================
# Rhombic Amplifier Examples (Section 3.1)
# ============================================================================

fn example_rhombic_single_point():
    """
    Single-point analysis of rhombic amplifier at theta=10deg.
    
    Reproduces the initial verification case from main_modular.mojo.
    Expected: ratio ~= cot(10deg) x 2 ~= 11.3 (single-sided)
    """
    print("\n" + "=" * 60)
    print("Rhombic Amplifier Single-Point Analysis (Section 3.1)")
    print("=" * 60)
    
    # Parameters matching paper Section 3.1
    var theta_deg = 10.0
    var geom = rhombic_default(
        L=15.0e-3, # 15mm total limb length
        h=1.0e-3, # 1mm flexure thickness
        d=10.0e-3, # 10mm depth
        theta_deg=theta_deg,
        L_flex=3.0e-3, # 3mm flexure hinges
    )
    
    var output_mass = _get_output_mass()
    var mat = aluminum_6061()
    var results = compute_results(
        mat, geom, output_mass,
        input_mass=0.0,
        single_sided=True
    )
    
    print("Parameters: theta =", theta_deg, "deg")
    print("Expected ratio (approx 2*cot theta):", 2.0 / tan(theta_deg * pi / 180.0))
    results.dump()


fn example_rhombic_sweep():
    """
    Theta sweep for rhombic amplifier.
    
    Reproduces Figure 4 from the paper with comparison curves
    from Refs [23], [28], and [31].
    
    Output: fig4_rhombic.png
    """
    print("\n" + "=" * 60)
    print("Generating Figure 4 (Rhombic Amplifier Sweep)...")
    print("=" * 60)
    
    # Sweep configuration matching main_modular.mojo
    var config = make_rhombic_config(
        L=15.0e-3,
        h=1.0e-3,
        d=10.0e-3,
        L_flex=3.0e-3,
    )
    
    var output_mass = _get_output_mass()
    var mat = aluminum_6061()
    
    # Run sweep
    var results = sweep_rhombic_theta(
        mat=mat,
        config=config,
        output_mass=output_mass,
        theta_start=0.01,
        theta_end=30.0,
        num_steps=500,
        single_sided=True,
        input_mass=0.0,
        verbose=True,
    )
    
    # Generate comparison curves
    var theta_values = get_theta_values(results)
    var comparisons = generate_rhombic_comparisons(
        mat,
        L=15.0e-3, # Use total limb length for comparison models
        h=1.0e-3,
        d=10.0e-3,
        theta_values=theta_values,
        single_sided=True,
    )
    
    # Plot with comparisons
    try:
        plot_rhombic_with_comparisons(
            results,
            comparisons,
            filename="fig4_rhombic.png",
        )
    except e:
        print("Plotting error:", e)


# ============================================================================
# Parallel Amplifier Examples (Section 3.2)
# ============================================================================

fn example_parallel_single_point():
    """
    Single-point analysis of parallel amplifier at equivalent theta=10deg.
    
    For parallel amplifier, H is derived from theta via H = L2 x tan(theta).
    """
    print("\n" + "=" * 60)
    print("Parallel Amplifier Single-Point Analysis (Section 3.2)")
    print("=" * 60)
    
    # Parameters matching paper Section 3.2
    var L_flex = 3.0e-3
    var L_link = 9.0e-3
    var h_flex = 0.5e-3
    var h_link = 4.0e-3
    var d = 10.0e-3
    
    # Compute H for equivalent theta = 10 deg
    var theta_target = 10.0
    var theta_rad = theta_target * pi / 180.0
    var H = L_link * tan(theta_rad)
    
    var geom = parallel_default(
        L_flex=L_flex,
        L_link=L_link,
        h_flex=h_flex,
        h_link=h_link,
        d=d,
        H=H,
    )
    
    var output_mass = _get_output_mass()
    var mat = aluminum_6061()
    var results = compute_results(
        mat, geom, output_mass,
        single_sided=True
    )
    
    print("Parameters: equivalent theta =", theta_target, "deg, H =", H * 1e3, "mm")
    results.dump()


fn example_parallel_sweep():
    """
    H sweep for parallel amplifier.
    
    Reproduces Figure 5 from the paper with comparison curves
    from Refs [23], [25], [26], and [31].
    
    Output: fig5_parallel.png
    """
    print("\n" + "=" * 60)
    print("Generating Figure 5 (Parallel Amplifier Sweep)...")
    print("=" * 60)
    
    # Sweep configuration matching main_modular.mojo
    var L_flex = 3.0e-3
    var L_link = 9.0e-3
    var h_flex = 0.5e-3
    var h_link = 4.0e-3
    var d = 10.0e-3
    
    var config = make_parallel_config(
        L_flex=L_flex,
        L_link=L_link,
        h_flex=h_flex,
        h_link=h_link,
        d=d,
    )
    
    var output_mass = _get_output_mass()
    var mat = aluminum_6061()
    
    # Run sweep
    var results = sweep_parallel_H(
        mat=mat,
        config=config,
        output_mass=output_mass,
        H_start=0.01e-3,
        H_end=5.0e-3,
        num_steps=500,
        single_sided=True,
        input_mass=0.0,
        verbose=True,
    )
    
    # Generate comparison curves
    var H_values = get_H_values_meters(results)
    var L0_total = L_flex + L_link # Total horizontal length for ref models
    
    var comparisons = generate_parallel_comparisons(
        mat,
        L0=L0_total,
        h_flex=h_flex,
        h_link=h_link,
        d=d,
        L_link=L_link,
        H_values=H_values,
        single_sided=True,
    )
    
    # Plot with comparisons
    try:
        plot_parallel_with_comparisons(
            results,
            comparisons,
            filename="fig5_parallel.png",
        )
    except e:
        print("Plotting error:", e)


# ============================================================================
# Aligned Amplifier Examples (Section 3.3)
# ============================================================================

fn example_aligned_sweep():
    """
    Theta sweep for aligned amplifier.
    
    Reproduces Figure 6 from the paper.
    
    Output: fig6_aligned.png
    """
    print("\n" + "=" * 60)
    print("Generating Figure 6 (Aligned Amplifier Sweep)...")
    print("=" * 60)
    
    # Sweep configuration matching main_modular.mojo
    var config = make_aligned_config(
        L_flex=3.0e-3,
        L_link=9.0e-3,
        h_flex=0.5e-3,
        h_link=4.0e-3,
        d=10.0e-3,
    )
    
    var output_mass = _get_output_mass()
    var mat = aluminum_6061()
    
    # Run sweep (fewer steps as in original)
    var results = sweep_aligned_theta(
        mat=mat,
        config=config,
        output_mass=output_mass,
        theta_start=0.01,
        theta_end=30.0,
        num_steps=100,
        single_sided=True,
        input_mass=0.0,
        verbose=True,
    )
    
    # Plot (no comparison models for aligned in the paper)
    try:
        plot_sweep_1d(results, "fig6_aligned.png", "Aligned Amplifier")
    except e:
        print("Plotting error:", e)


# ============================================================================
# Run All Examples
# ============================================================================

fn run_all_examples():
    """
    Run all example analyses.
    
    Produces:
    - Console output with single-point results
    - fig4_rhombic.png
    - fig5_parallel.png
    - fig6_aligned.png
    
    These can be compared against the original main_modular.mojo outputs.
    """
    print("\n" + "#" * 70)
    print("#  FLEXURE AMPLIFIER LIBRARY - VALIDATION EXAMPLES")
    print("#" * 70)
    
    # Single-point analyses
    example_rhombic_single_point()
    example_parallel_single_point()
    
    # Sweep analyses (produce plots)
    example_rhombic_sweep()
    example_parallel_sweep()
    example_aligned_sweep()
    
    print("\n" + "#" * 70)
    print("#  ALL EXAMPLES COMPLETE")
    print("#" * 70)
    print("\nGenerated plots:")
    print("  - fig4_rhombic.png  (Rhombic amplifier, Figure 4)")
    print("  - fig5_parallel.png (Parallel amplifier, Figure 5)")
    print("  - fig6_aligned.png  (Aligned amplifier, Figure 6)")
