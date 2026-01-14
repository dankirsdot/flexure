"""
Piezo Amplifier Design Exploration Tool

Explores custom amplifier designs for a single SCL070718 piezo stack:
- Stack: 7mm x 7mm x 18mm, stroke ~22um, blocked force ~1900N
- Configuration: Single stack centrally between input ports
- Negative theta for preload and correct motion direction

Physical rationale:
- Negative theta (limbs lean inward) provides compressive preload on the piezo
- Piezo expansion pushes inputs apart -> output extends upward (positive motion)
- single_sided=True: Both inputs move symmetrically from central actuator

Performs 2D sweeps over theta and flexure thickness for all three amplifier
types and multiple materials, generating contour plots and a ranked summary.

Usage:
    mojo run design.mojo
"""

from math import pi, sin, cos, tan, sqrt, atan2
from python import Python, PythonObject

from model import MaterialProperties, GeometricParameters
from materials import (
    aluminum_6061, aluminum_7075,
    steel_304ss, steel_17_4ph,
    titanium_ti6al4v,
    beryllium_copper, invar_36
)
from geometry import rhombic, parallel, aligned, InputPortParams
from analysis import compute_results, make_output_mass_from_block, OutputMass


# ============================================================================
# Piezo Stack Parameters (SCL070718 from PiezoDrive)
# ============================================================================

comptime PIEZO_LENGTH: Float64 = 18.0e-3 # m - stack length (actuation axis)
comptime PIEZO_WIDTH: Float64 = 7.0e-3 # m - stack cross-section
comptime PIEZO_STROKE: Float64 = 22.0e-6 # m - nominal bipolar stroke
comptime PIEZO_BLOCKED_FORCE: Float64 = 1900.0 # N
comptime PIEZO_STIFFNESS: Float64 = 100.0e6 # N/m (~100 N/um)


# ============================================================================
# Base Geometry Configuration (piezo-scaled)
# ============================================================================

# Input port: Accommodates piezo stack with clearance
comptime L0_BASE: Float64 = 12.0e-3 # m - vertical input port height
comptime H0_BASE: Float64 = 5.0e-3 # m - horizontal input port width
comptime D_BASE: Float64 = 9.0e-3 # m - depth (matches stack width + margin)

# D_BASE sweep parameters
comptime D_BASE_MIN: Float64 = 8.0e-3 # m (8 mm)
comptime D_BASE_MAX: Float64 = 10.0e-3 # m (12 mm)
comptime D_BASE_STEP: Float64 = 0.5e-3 # m (0.5 mm step)

# Limb dimensions
comptime L_TOTAL: Float64 = 15.0e-3 # m - total limb length
comptime L_FLEX: Float64 = 3.0e-3 # m - flexure hinge length
comptime L_LINK: Float64 = 9.0e-3 # m - rigid link length (for parallel/aligned)
comptime H_LINK: Float64 = 4.0e-3 # m - link thickness (essentially rigid)


# ============================================================================
# Sweep Configuration
# ============================================================================

comptime THETA_MIN: Float64 = -30.0 # degrees (most negative)
comptime THETA_MAX: Float64 = -0.5 # degrees (avoid singularity near 0)
comptime THETA_STEPS: Int = 100

comptime H_FLEX_MIN: Float64 = 1.0e-3 # m (1 mm)
comptime H_FLEX_MAX: Float64 = 5.0e-3 # m (5 mm)
comptime H_FLEX_STEPS: Int = 100


# ============================================================================
# Design Constraints
# ============================================================================

comptime RATIO_MIN: Float64 = 10.0 # Minimum amplification ratio
comptime STIFFNESS_MAX: Float64 = 50.0 # Maximum input stiffness (N/um)
comptime PRELOAD_MIN: Float64 = 200.0 # Minimum preload force (N)


# ============================================================================
# Helper: Create output mass (small platform)
# ============================================================================

fn get_design_output_mass() -> OutputMass:
    """Output mass: Small 5x5x8 mm aluminum block (~0.5g)."""
    return make_output_mass_from_block(5.0e-3, 5.0e-3, 8.0e-3, 2770.0)


fn get_input_port() -> InputPortParams:
    """Input port sized for piezo stack."""
    return InputPortParams(L0=L0_BASE, h0=H0_BASE)


# ============================================================================
# Geometry Factories for Negative Theta
# ============================================================================

fn make_rhombic_geometry(theta_deg: Float64, h_flex: Float64, d: Float64) -> GeometricParameters:
    """
    Create rhombic geometry with negative theta.
    
    Uses uniform thickness for flexures and link (classic rhombic).
    """
    return rhombic(
        L=L_TOTAL,
        h=h_flex,
        d=d,
        theta_deg=theta_deg,
        L_flex=L_FLEX,
        input_port=get_input_port(),
    )


fn make_parallel_geometry(theta_deg: Float64, h_flex: Float64, d: Float64) -> GeometricParameters:
    """
    Create parallel geometry with negative theta.
    
    For parallel type, H (vertical offset) derives from theta:
    H = L_link x tan(|theta|), but negative theta means H is conceptually
    "below" the horizontal. The geometry factory handles signs internally.
    """
    var theta_rad = theta_deg * pi / 180.0
    var H = L_LINK * tan(-theta_rad) # Positive H for negative theta
    
    return parallel(
        L_flex=L_FLEX,
        L_link=L_LINK,
        h_flex=h_flex,
        h_link=H_LINK,
        d=d,
        H=H,
        input_port=get_input_port(),
    )


fn make_aligned_geometry(theta_deg: Float64, h_flex: Float64, d: Float64) -> GeometricParameters:
    """
    Create aligned geometry with negative theta.
    
    All beams share the same inclination angle.
    """
    return aligned(
        L_flex=L_FLEX,
        L_link=L_LINK,
        h_flex=h_flex,
        h_link=H_LINK,
        d=d,
        theta_deg=theta_deg,
        input_port=get_input_port(),
    )


# ============================================================================
# Preload Estimation
# ============================================================================

fn estimate_preload(theta_deg: Float64, k_in: Float64) -> Float64:
    """
    Estimate piezo preload force from negative-theta geometry.
    
    When theta < 0, the V-shaped limbs want to squeeze inward.
    The preload displacement is approximately:
        delta ~= 2 x L_total x |sin(theta)|
    
    And preload force:
        F_preload ~= k_in x delta
    
    This is a rough approximation - actual preload depends on assembly.
    
    Args:
        theta_deg: Inclination angle (should be negative).
        k_in: Input stiffness in N/m.
    
    Returns:
        Estimated preload force in N.
    """
    var theta_rad = theta_deg * pi / 180.0
    var delta = 2.0 * L_TOTAL * sin(-theta_rad) # Positive for negative theta
    return k_in * delta


# ============================================================================
# 2D Sweep for Single Structure/Material
# ============================================================================

@fieldwise_init
struct DesignResult(Copyable, Movable):
    """Result from a single design point."""
    var theta_deg: Float64
    var h_flex_mm: Float64
    var ratio: Float64
    var stiffness: Float64 # N/um
    var f1: Float64
    var f2: Float64
    var f3: Float64
    var preload: Float64 # N
    var meets_constraints: Bool
    
    fn __init__(out self):
        self.theta_deg = 0.0
        self.h_flex_mm = 0.0
        self.ratio = 0.0
        self.stiffness = 0.0
        self.f1 = 0.0
        self.f2 = 0.0
        self.f3 = 0.0
        self.preload = 0.0
        self.meets_constraints = False
    
    fn __copyinit__(out self, existing: Self):
        self.theta_deg = existing.theta_deg
        self.h_flex_mm = existing.h_flex_mm
        self.ratio = existing.ratio
        self.stiffness = existing.stiffness
        self.f1 = existing.f1
        self.f2 = existing.f2
        self.f3 = existing.f3
        self.preload = existing.preload
        self.meets_constraints = existing.meets_constraints
    
    fn __moveinit__(out self, deinit existing: Self):
        self.theta_deg = existing.theta_deg
        self.h_flex_mm = existing.h_flex_mm
        self.ratio = existing.ratio
        self.stiffness = existing.stiffness
        self.f1 = existing.f1
        self.f2 = existing.f2
        self.f3 = existing.f3
        self.preload = existing.preload
        self.meets_constraints = existing.meets_constraints
    
    fn copy(self) -> Self:
        """Create an explicit copy."""
        var result = DesignResult()
        result.theta_deg = self.theta_deg
        result.h_flex_mm = self.h_flex_mm
        result.ratio = self.ratio
        result.stiffness = self.stiffness
        result.f1 = self.f1
        result.f2 = self.f2
        result.f3 = self.f3
        result.preload = self.preload
        result.meets_constraints = self.meets_constraints
        return result^


fn run_2d_sweep(
    mat: MaterialProperties,
    structure: String,
    d: Float64,
    verbose: Bool = True,
) -> List[DesignResult]:
    """
    Run 2D parameter sweep over theta and h_flex.
    
    Args:
        mat: Material properties.
        structure: "rhombic", "parallel", or "aligned".
        d: Depth dimension in meters.
        verbose: Print progress updates.
    
    Returns:
        List of DesignResult for all grid points.
    """
    var results = List[DesignResult]()
    var output_mass = get_design_output_mass()
    
    var theta_step = (THETA_MAX - THETA_MIN) / Float64(THETA_STEPS - 1)
    var h_step = (H_FLEX_MAX - H_FLEX_MIN) / Float64(H_FLEX_STEPS - 1)
    
    for i in range(THETA_STEPS):
        var theta_deg = THETA_MIN + Float64(i) * theta_step
        
        for j in range(H_FLEX_STEPS):
            var h_flex = H_FLEX_MIN + Float64(j) * h_step
            
            # Create geometry based on structure type
            var geom: GeometricParameters
            if structure == "rhombic":
                geom = make_rhombic_geometry(theta_deg, h_flex, d)
            elif structure == "parallel":
                geom = make_parallel_geometry(theta_deg, h_flex, d)
            else: # aligned
                geom = make_aligned_geometry(theta_deg, h_flex, d)
            
            # Compute analysis
            var analysis = compute_results(
                mat, geom, output_mass,
                input_mass=0.0,
                single_sided=True
            )
            
            # Store result
            var result = DesignResult()
            result.theta_deg = theta_deg
            result.h_flex_mm = h_flex * 1000.0 # Convert to mm
            result.ratio = analysis.amplification_ratio
            result.stiffness = analysis.input_stiffness / 1.0e6 # N/um
            result.f1 = analysis.f1()
            result.f2 = analysis.f2()
            result.f3 = analysis.f3()
            result.preload = estimate_preload(theta_deg, analysis.input_stiffness)
            
            # Check constraints
            result.meets_constraints = (
                result.ratio >= RATIO_MIN and
                result.stiffness <= STIFFNESS_MAX and
                result.preload >= PRELOAD_MIN
            )
            
            results.append(result^)
        
        if verbose and i % 10 == 0:
            print("  Progress:", i + 1, "/", THETA_STEPS, "theta values")
    
    return results^


# ============================================================================
# Plotting Functions
# ============================================================================

fn plot_ratio_stiffness_contours(
    results: List[DesignResult],
    structure: String,
    material_name: String,
    d_mm: Float64,
) raises:
    """
    Generate 2-panel contour plot for ratio and stiffness.
    
    Panels: (a) Amplification Ratio, (b) Input Stiffness
    Both plots are square aspect ratio.
    """
    var plt = Python.import_module("matplotlib.pyplot")
    var np = Python.import_module("numpy")
    var builtins = Python.import_module("builtins")
    
    # Create numpy arrays for contour plot
    var theta_vals = np.linspace(THETA_MIN, THETA_MAX, THETA_STEPS)
    var h_vals = np.linspace(H_FLEX_MIN * 1000.0, H_FLEX_MAX * 1000.0, H_FLEX_STEPS)
    
    # Build 2D data as nested Python lists
    var ratio_rows = builtins.list()
    var stiff_rows = builtins.list()
    
    for i in range(THETA_STEPS):
        var ratio_row = builtins.list()
        var stiff_row = builtins.list()
        
        for j in range(H_FLEX_STEPS):
            var idx = i * H_FLEX_STEPS + j
            ratio_row.append(results[idx].ratio)
            stiff_row.append(results[idx].stiffness)
        
        ratio_rows.append(ratio_row)
        stiff_rows.append(stiff_row)
    
    # Convert to numpy arrays
    var ratio_grid = np.abs(np.array(ratio_rows))
    var stiff_grid = np.array(stiff_rows)
    
    # Create figure with 2 square subplots (1 row x 2 columns)
    var fig_info = plt.subplots(1, 2)
    var fig = fig_info[0]
    var axes = fig_info[1]
    fig.set_size_inches(10, 5)
    
    var d_str = String(d_mm)
    fig.suptitle(structure.upper() + " Amplifier - " + material_name + " (d=" + d_str + "mm)", fontsize=14)
    
    # Subplot (a): Amplification Ratio
    var ax1 = axes[0]
    var c1 = ax1.contourf(theta_vals, h_vals, ratio_grid.T, 20)
    plt.colorbar(c1, ax=ax1)
    ax1.set_xlabel(String("theta (deg)"))
    ax1.set_ylabel(String("h_flex (mm)"))
    ax1.set_title(String("(a) Amplification Ratio"))
    
    # Subplot (b): Input Stiffness
    var ax2 = axes[1]
    var c2 = ax2.contourf(theta_vals, h_vals, stiff_grid.T, 20)
    plt.colorbar(c2, ax=ax2)
    ax2.set_xlabel(String("theta (deg)"))
    ax2.set_ylabel(String("h_flex (mm)"))
    ax2.set_title(String("(b) Input Stiffness (N/um)"))
    
    plt.tight_layout()
    
    var filename = structure + "_" + material_name.replace(" ", "_") + "_d_" + d_str + "mm_mechanics.png"
    plt.savefig(filename, dpi=150)
    print("Saved:", filename)
    _ = plt.close(fig)


fn plot_frequency_contours(
    results: List[DesignResult],
    structure: String,
    material_name: String,
    d_mm: Float64,
) raises:
    """
    Generate 3-panel contour plot for natural frequencies.
    
    Panels: (a) f1, (b) f2, (c) f3
    Horizontal layout (1 row x 3 columns).
    """
    var plt = Python.import_module("matplotlib.pyplot")
    var np = Python.import_module("numpy")
    var builtins = Python.import_module("builtins")
    
    # Create numpy arrays for contour plot
    var theta_vals = np.linspace(THETA_MIN, THETA_MAX, THETA_STEPS)
    var h_vals = np.linspace(H_FLEX_MIN * 1000.0, H_FLEX_MAX * 1000.0, H_FLEX_STEPS)
    
    # Build 2D data as nested Python lists
    var f1_rows = builtins.list()
    var f2_rows = builtins.list()
    var f3_rows = builtins.list()
    
    for i in range(THETA_STEPS):
        var f1_row = builtins.list()
        var f2_row = builtins.list()
        var f3_row = builtins.list()
        
        for j in range(H_FLEX_STEPS):
            var idx = i * H_FLEX_STEPS + j
            f1_row.append(results[idx].f1)
            f2_row.append(results[idx].f2)
            f3_row.append(results[idx].f3)
        
        f1_rows.append(f1_row)
        f2_rows.append(f2_row)
        f3_rows.append(f3_row)
    
    # Convert to numpy arrays
    var f1_grid = np.array(f1_rows)
    var f2_grid = np.array(f2_rows)
    var f3_grid = np.array(f3_rows)
    
    # Create figure with 3 subplots (1 row x 3 columns)
    var fig_info = plt.subplots(1, 3)
    var fig = fig_info[0]
    var axes = fig_info[1]
    fig.set_size_inches(15, 5)
    
    var d_str = String(d_mm)
    fig.suptitle(structure.upper() + " Amplifier - " + material_name + " (d=" + d_str + "mm) - Natural Frequencies", fontsize=14)
    
    # Subplot (a): f1
    var ax1 = axes[0]
    var c1 = ax1.contourf(theta_vals, h_vals, f1_grid.T, 20)
    plt.colorbar(c1, ax=ax1)
    ax1.set_xlabel(String("theta (deg)"))
    ax1.set_ylabel(String("h_flex (mm)"))
    ax1.set_title(String("(a) f1 (Hz)"))
    
    # Subplot (b): f2
    var ax2 = axes[1]
    var c2 = ax2.contourf(theta_vals, h_vals, f2_grid.T, 20)
    plt.colorbar(c2, ax=ax2)
    ax2.set_xlabel(String("theta (deg)"))
    ax2.set_ylabel(String("h_flex (mm)"))
    ax2.set_title(String("(b) f2 (Hz)"))
    
    # Subplot (c): f3
    var ax3 = axes[2]
    var c3 = ax3.contourf(theta_vals, h_vals, f3_grid.T, 20)
    plt.colorbar(c3, ax=ax3)
    ax3.set_xlabel(String("theta (deg)"))
    ax3.set_ylabel(String("h_flex (mm)"))
    ax3.set_title(String("(c) f3 (Hz)"))
    
    plt.tight_layout()
    
    var filename = structure + "_" + material_name.replace(" ", "_") + "_d_" + d_str + "mm_frequencies.png"
    plt.savefig(filename, dpi=150)
    print("Saved:", filename)
    _ = plt.close(fig)


# ============================================================================
# Ranking and Summary
# ============================================================================

fn find_best_configurations(
    all_results: List[DesignResult],
    structure: String,
    material: String,
    top_n: Int = 5,
) -> List[DesignResult]:
    """
    Find top N configurations that meet constraints, sorted by ratio.
    """
    var valid = List[DesignResult]()
    
    for i in range(len(all_results)):
        if all_results[i].meets_constraints:
            valid.append(all_results[i].copy())
    
    # Simple bubble sort by ratio (descending) - Mojo lacks sorting
    for i in range(len(valid)):
        for j in range(i + 1, len(valid)):
            if valid[j].ratio > valid[i].ratio:
                var temp = valid[i].copy()
                valid[i] = valid[j].copy()
                valid[j] = temp^
    
    # Return top N
    var result = List[DesignResult]()
    for i in range(min(top_n, len(valid))):
        result.append(valid[i].copy())
    
    return result^


fn print_summary_table(
    rankings: List[DesignResult],
    structures: List[String],
    materials: List[String],
):
    """Print ranked summary table."""
    print("\n" + "=" * 90)
    print("TOP CONFIGURATIONS (meeting constraints: ratio >=", RATIO_MIN, 
          ", k_in <=", STIFFNESS_MAX, "N/um, preload >=", PRELOAD_MIN, "N)")
    print("=" * 90)
    print("Rank | Structure | Material    | theta(deg) | h(mm) | Ratio | k_in  | f1(Hz) | Preload(N)")
    print("-" * 90)
    
    for i in range(len(rankings)):
        var r = rankings[i].copy()
        print(
            "  ", i + 1, " |",
            "          |", # Will need to add structure/material tracking
            "            |",
            "   ", r.theta_deg, " |",
            " ", r.h_flex_mm, " |",
            " ", r.ratio, " |",
            " ", r.stiffness, " |",
            " ", r.f1, " |",
            " ", r.preload
        )
    
    print("=" * 90)


# ============================================================================
# Main Entry Point
# ============================================================================

fn main():
    """Run full design exploration."""
    print("\n" + "#" * 70)
    print("#  PIEZO AMPLIFIER DESIGN EXPLORATION")
    print("#  Stack: SCL070718 (7x7x18mm, ~24um stroke)")
    print("#" * 70)
    
    # Define materials and structures to explore
    var structures = List[String]()
    # structures.append("rhombic")
    # structures.append("parallel")
    structures.append("aligned")
    
    # Material info (name, properties pairs)
    var mat_names = List[String]()
    mat_names.append("Al6061")
    # mat_names.append("Al7075")
    # mat_names.append("SS304")
    # mat_names.append("SS17-4PH")
    # mat_names.append("Ti6Al4V")
    # mat_names.append("BeCu")
    # mat_names.append("Invar36")
    
    # var output_mass = get_design_output_mass()
    
    # Collect all top results for final ranking
    var all_top_results = List[DesignResult]()
    var total_plots = 0
    
    # Outer loop over D_BASE values
    var d = D_BASE_MIN
    while d <= D_BASE_MAX + 1e-6:
        var d_mm = d * 1000.0
        print("\n" + "=" * 70)
        print("D_BASE =", d_mm, "mm")
        print("=" * 70)
        
        # Run sweeps for each combination
        for s_idx in range(len(structures)):
            var structure = structures[s_idx]
            
            for m_idx in range(len(mat_names)):
                var mat_name = mat_names[m_idx]
                var mat: MaterialProperties
                
                if mat_name == "Al6061":
                    mat = aluminum_6061()
                elif mat_name == "Al7075":
                    mat = aluminum_7075()
                elif mat_name == "SS304":
                    mat = steel_304ss()
                elif mat_name == "SS17-4PH":
                    mat = steel_17_4ph()
                elif mat_name == "Ti6Al4V":
                    mat = titanium_ti6al4v()
                elif mat_name == "BeCu":
                    mat = beryllium_copper()
                else:
                    mat = invar_36()
                
                print("\n" + "-" * 50)
                print("Running:", structure, "-", mat_name, "- d =", d_mm, "mm")
                print("-" * 50)
                
                # Run 2D sweep
                var results = run_2d_sweep(mat, structure, d, verbose=True)
                
                # Generate contour plots (2 separate files)
                try:
                    plot_ratio_stiffness_contours(results, structure, mat_name, d_mm)
                    plot_frequency_contours(results, structure, mat_name, d_mm)
                    total_plots += 2
                except e:
                    print("Plotting error:", e)
                
                # Find best configurations
                var best = find_best_configurations(results, structure, mat_name, top_n=3)
                print("Found", len(best), "configurations meeting constraints")
                
                # Add to global ranking
                for i in range(len(best)):
                    all_top_results.append(best[i].copy())
        
        d += D_BASE_STEP
    
    # Print summary
    print("\n" + "#" * 70)
    print("#  DESIGN EXPLORATION COMPLETE")
    print("#" * 70)
    print("\nGenerated", total_plots, "contour plot files")
    print("Total configurations meeting constraints:", len(all_top_results))
