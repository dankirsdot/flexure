"""
Centralized plotting functions for flexure amplifier analysis.

This module provides reusable plotting functions that work with sweep results
from any amplifier type. It eliminates the duplicated plotting code from
the original main_modular.mojo.

All plotting uses matplotlib through Python interop.

Usage:
    from plotting import plot_sweep_1d, plot_rhombic_with_comparisons
    
    let results = sweep_rhombic_theta(...)
    plot_sweep_1d(results, "fig4_rhombic.png", title_prefix="Rhombic Amplifier")
"""

from python import Python
from sweeps import Sweep1DResults
from references import RhombicComparisonResults, ParallelComparisonResults


# ============================================================================
# Core Plotting Function
# ============================================================================

fn plot_sweep_1d(
    results: Sweep1DResults,
    filename: String,
    title_prefix: String = "",
) raises:
    """
    Generate the standard 3-subplot figure (ratio, stiffness, frequencies).
    
    This is the basic plotting function that works with any 1D sweep results.
    For plots with comparison curves, use the type-specific functions below.
    
    Args:
        results: Sweep1DResults from any sweep function.
        filename: Output filename (e.g., "fig6_aligned.png").
        title_prefix: Optional prefix for subplot titles.
    
    Raises:
        If matplotlib is not available or plotting fails.
    """
    var plt = Python.import_module("matplotlib.pyplot")
    var builtins = Python.import_module("builtins")
    
    # Convert Mojo lists to Python lists
    var x_vals = builtins.list()
    var ratios = builtins.list()
    var stiffnesses = builtins.list()
    var f1_vals = builtins.list()
    var f2_vals = builtins.list()
    var f3_vals = builtins.list()
    
    for i in range(len(results.param_values)):
        x_vals.append(results.param_values[i])
        ratios.append(results.ratios[i])
        stiffnesses.append(results.stiffnesses[i])
        f1_vals.append(results.f1[i])
        f2_vals.append(results.f2[i])
        f3_vals.append(results.f3[i])
    
    # Create figure
    var fig = plt.figure()
    fig.set_size_inches(8, 12)
    
    # Build axis label
    var x_label = results.param_name + " (" + results.param_unit + ")"
    
    # Subplot (a): Amplification Ratio
    var ax1 = plt.subplot(3, 1, 1)
    ax1.plot(x_vals, ratios, "r-")
    ax1.set_ylabel("Amplification Ratio")
    if title_prefix:
        ax1.set_title("(a) " + title_prefix + " - Displacement amplification ratio")
    else:
        ax1.set_title("(a) Displacement amplification ratio")
    ax1.grid(True)
    
    var legend_a = builtins.list()
    legend_a.append("Model")
    ax1.legend(legend_a)
    
    # Subplot (b): Input Stiffness
    var ax2 = plt.subplot(3, 1, 2)
    ax2.plot(x_vals, stiffnesses, "r-")
    ax2.set_ylabel("Input Stiffness (N/um)")
    if title_prefix:
        ax2.set_title("(b) " + title_prefix + " - Input stiffness")
    else:
        ax2.set_title("(b) Input stiffness")
    ax2.grid(True)
    
    # Subplot (c): Natural Frequencies
    var ax3 = plt.subplot(3, 1, 3)
    ax3.plot(x_vals, f1_vals, "r-")
    ax3.plot(x_vals, f2_vals, "g--")
    ax3.plot(x_vals, f3_vals, "b-.")
    ax3.set_ylabel("Natural Frequency (Hz)")
    ax3.set_xlabel(x_label)
    if title_prefix:
        ax3.set_title("(c) " + title_prefix + " - Natural frequencies")
    else:
        ax3.set_title("(c) The first three in-plane natural frequencies")
    ax3.grid(True)
    
    var legend_c = builtins.list()
    legend_c.append("1st Mode")
    legend_c.append("2nd Mode")
    legend_c.append("3rd Mode")
    ax3.legend(legend_c)
    
    plt.tight_layout()
    plt.savefig(filename)
    print("Plot saved to", filename)
    _ = plt.close(fig)


# ============================================================================
# Rhombic Plotting with Comparison Curves
# ============================================================================

fn plot_rhombic_with_comparisons(
    results: Sweep1DResults,
    comparisons: RhombicComparisonResults,
    filename: String,
) raises:
    """
    Plot rhombic amplifier sweep with literature comparison curves.
    
    Reproduces Figure 4 from the paper with the proposed model and
    three reference models (Refs 23, 28, 31).
    
    Args:
        results: Sweep1DResults from sweep_rhombic_theta.
        comparisons: RhombicComparisonResults from generate_rhombic_comparisons.
        filename: Output filename.
    """
    var plt = Python.import_module("matplotlib.pyplot")
    var builtins = Python.import_module("builtins")
    
    # Convert Mojo lists to Python lists
    var x_vals = builtins.list()
    var ratios = builtins.list()
    var stiffnesses = builtins.list()
    var f1_vals = builtins.list()
    var f2_vals = builtins.list()
    var f3_vals = builtins.list()
    var r23_py = builtins.list()
    var r28_py = builtins.list()
    var r31_py = builtins.list()
    
    for i in range(len(results.param_values)):
        x_vals.append(results.param_values[i])
        ratios.append(results.ratios[i])
        stiffnesses.append(results.stiffnesses[i])
        f1_vals.append(results.f1[i])
        f2_vals.append(results.f2[i])
        f3_vals.append(results.f3[i])
        r23_py.append(comparisons.ref23[i])
        r28_py.append(comparisons.ref28[i])
        r31_py.append(comparisons.ref31[i])
    
    # Create figure
    var fig = plt.figure()
    fig.set_size_inches(8, 12)
    
    # Subplot (a): Ratio with comparisons
    var ax1 = plt.subplot(3, 1, 1)
    ax1.plot(x_vals, ratios, "r-")
    ax1.plot(x_vals, r23_py, "k-")
    ax1.plot(x_vals, r28_py, "g--")
    ax1.plot(x_vals, r31_py, "b-.")
    ax1.set_ylabel("Amplification Ratio")
    ax1.set_title("(a) Displacement amplification ratio")
    ax1.grid(True)
    
    var legend_a = builtins.list()
    legend_a.append("Proposed Model")
    legend_a.append("Ref [23]")
    legend_a.append("Ref [28]")
    legend_a.append("Ref [31]")
    ax1.legend(legend_a)
    
    # Subplot (b): Input Stiffness
    var ax2 = plt.subplot(3, 1, 2)
    ax2.plot(x_vals, stiffnesses, "r-")
    ax2.set_ylabel("Input Stiffness (N/um)")
    ax2.set_title("(b) Input stiffness")
    ax2.grid(True)
    
    # Subplot (c): Natural Frequencies
    var ax3 = plt.subplot(3, 1, 3)
    ax3.plot(x_vals, f1_vals, "r-")
    ax3.plot(x_vals, f2_vals, "g--")
    ax3.plot(x_vals, f3_vals, "b-.")
    ax3.set_ylabel("Natural Frequency (Hz)")
    ax3.set_xlabel("Angle (deg)")
    ax3.set_title("(c) The first three in-plane natural frequencies")
    ax3.grid(True)
    
    var legend_c = builtins.list()
    legend_c.append("1st Mode")
    legend_c.append("2nd Mode")
    legend_c.append("3rd Mode")
    ax3.legend(legend_c)
    
    plt.tight_layout()
    plt.savefig(filename)
    print("Plot saved to", filename)
    _ = plt.close(fig)


# ============================================================================
# Parallel Plotting with Comparison Curves
# ============================================================================

fn plot_parallel_with_comparisons(
    results: Sweep1DResults,
    comparisons: ParallelComparisonResults,
    filename: String,
) raises:
    """
    Plot parallel amplifier sweep with literature comparison curves.
    
    Reproduces Figure 5 from the paper with the proposed model and
    four reference models (Refs 23, 25, 26, 31).
    
    Args:
        results: Sweep1DResults from sweep_parallel_H.
        comparisons: ParallelComparisonResults from generate_parallel_comparisons.
        filename: Output filename.
    """
    var plt = Python.import_module("matplotlib.pyplot")
    var builtins = Python.import_module("builtins")
    
    # Convert to Python lists
    var x_vals = builtins.list()
    var ratios = builtins.list()
    var stiffnesses = builtins.list()
    var f1_vals = builtins.list()
    var f2_vals = builtins.list()
    var f3_vals = builtins.list()
    var r23_py = builtins.list()
    var r25_py = builtins.list()
    var r26_py = builtins.list()
    var r31_py = builtins.list()
    
    for i in range(len(results.param_values)):
        x_vals.append(results.param_values[i])
        ratios.append(results.ratios[i])
        stiffnesses.append(results.stiffnesses[i])
        f1_vals.append(results.f1[i])
        f2_vals.append(results.f2[i])
        f3_vals.append(results.f3[i])
        r23_py.append(comparisons.ref23[i])
        r25_py.append(comparisons.ref25[i])
        r26_py.append(comparisons.ref26[i])
        r31_py.append(comparisons.ref31[i])
    
    # Create figure
    var fig = plt.figure()
    fig.set_size_inches(8, 12)
    
    # Subplot (a): Ratio with comparisons
    var ax1 = plt.subplot(3, 1, 1)
    ax1.plot(x_vals, ratios, "r-")
    ax1.plot(x_vals, r23_py, "k-")
    ax1.plot(x_vals, r25_py, "c--")
    ax1.plot(x_vals, r26_py, "b-.")
    ax1.plot(x_vals, r31_py, "m:")
    ax1.set_ylabel("Disp Ratio")
    ax1.set_title("(a) Displacement amplification ratio")
    ax1.grid(True)
    
    var legend_a = builtins.list()
    legend_a.append("Model")
    legend_a.append("Ref 23")
    legend_a.append("Ref 25")
    legend_a.append("Ref 26")
    legend_a.append("Ref 31")
    ax1.legend(legend_a)
    
    # Subplot (b): Input Stiffness
    var ax2 = plt.subplot(3, 1, 2)
    ax2.plot(x_vals, stiffnesses, "r-")
    ax2.set_ylabel("Input Stiffness (N/um)")
    ax2.set_title("(b) Input stiffness")
    ax2.grid(True)
    
    # Subplot (c): Natural Frequencies
    var ax3 = plt.subplot(3, 1, 3)
    ax3.plot(x_vals, f1_vals, "r-")
    ax3.plot(x_vals, f2_vals, "g--")
    ax3.plot(x_vals, f3_vals, "b-.")
    ax3.set_ylabel("Natural Frequency (Hz)")
    ax3.set_xlabel("H (mm)")
    ax3.set_title("(c) Natural frequencies")
    ax3.grid(True)
    
    var legend_c = builtins.list()
    legend_c.append("1st Mode")
    legend_c.append("2nd Mode")
    legend_c.append("3rd Mode")
    ax3.legend(legend_c)
    
    plt.tight_layout()
    plt.savefig(filename)
    print("Plot saved to", filename)
    _ = plt.close(fig)
