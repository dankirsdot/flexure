from model import *
from python import Python
from math import pi, sin, cos, tan, sqrt, atan2

struct AmplifierResults(Copyable):
    var amplification_ratio: Float64
    var input_stiffness: Float64
    var natural_frequencies: List[Float64]
    
    fn __init__(out self):
        self.amplification_ratio = 0.0
        self.input_stiffness = 0.0
        self.natural_frequencies = List[Float64]()
        
    fn __copyinit__(out self, existing: Self):
        self.amplification_ratio = existing.amplification_ratio
        self.input_stiffness = existing.input_stiffness
        self.natural_frequencies = List[Float64](capacity=len(existing.natural_frequencies))
        for i in range(len(existing.natural_frequencies)):
            self.natural_frequencies.append(existing.natural_frequencies[i])

    fn dump(self):
        print("Amplification Ratio: ", self.amplification_ratio)
        print("Input Stiffness: ", self.input_stiffness)
        print("Natural Frequencies: ")
        for i in range(len(self.natural_frequencies)):
            print("  Mode ", i+1, ": ", self.natural_frequencies[i], " Hz")

fn get_rhombic_geometry(
    L: Float64, h: Float64,
    L0: Float64, h0: Float64,
    d: Float64,
    theta_deg: Float64,
    L_flex: Float64
) -> GeometricParameters:
    var theta_rad = theta_deg * pi / 180.0
    var L_rigid = L - 2.0 * L_flex
    
    # Beam 1 (Input Port Half)
    var L_1 = sqrt((L0/2)**2 + (h0/2)**2)
    var h_1 = h0 
    var theta_1 = atan2(L0/2, h0/2)

    # Beam 2 (Flexure 1)
    var L_2 = L_flex
    var h_2 = h
    var theta_2 = theta_rad

    # Beam 3 (Rigid Link)
    var L_3 = L_rigid
    var h_3 = h
    var theta_3 = theta_rad

    # Beam 4 (Flexure 2)
    var L_4 = L_flex
    var h_4 = h
    var theta_4 = theta_rad

    return GeometricParameters(
        L_1=L_1, L_2=L_2, L_3=L_3, L_4=L_4,
        h_1=h_1, h_2=h_2, h_3=h_3, h_4=h_4,
        theta_1=theta_1, theta_2=theta_2, theta_3=theta_3, theta_4=theta_4,
        H=0.0, d=d
    )

fn solve_rhombic_amplifier(
    mat: MaterialProperties,
    L: Float64, h: Float64,
    L0: Float64, h0: Float64,
    d: Float64,
    theta_deg: Float64,
    L_flex: Float64, 
    m_out: Float64,
    J_out: Float64,
    m_in: Float64,
    use_single_sided: Bool
) -> AmplifierResults:
    var geom = get_rhombic_geometry(L, h, L0, h0, d, theta_deg, L_flex)
    
    var res = AmplifierResults()
    res.amplification_ratio = calculate_amplification_ratio(mat, geom, use_single_sided)
    res.input_stiffness = calculate_input_stiffness(mat, geom, use_single_sided)
    res.natural_frequencies = calculate_natural_frequencies(mat, geom, m_out, J_out, m_in)
    
    return res^

fn run_rhombic_sweep(
    mat: MaterialProperties,
    L: Float64, h: Float64,
    L0: Float64, h0: Float64,
    d: Float64,
    L_flex: Float64,
    m_out: Float64, J_out: Float64, m_in: Float64,
    use_single_sided: Bool,
    num_steps: Int,
    theta_start: Float64,
    theta_end: Float64
):
    print("\n------------------------------------------------")
    print("Generating Figure 4 (a, b, c) Sweep...")
    
    try:
        var plt = Python.import_module("matplotlib.pyplot")
        var np = Python.import_module("numpy")
        var builtins = Python.import_module("builtins")
        
        var thetas_deg = builtins.list()
        var R_model = builtins.list()
        var R_23 = builtins.list()
        var R_28 = builtins.list()
        var R_31 = builtins.list()
        var K_in_model = builtins.list()
        var f1_model = builtins.list()
        var f2_model = builtins.list()
        var f3_model = builtins.list()
        
        # Theoretical Constants for Comparison Models
        var A = d * h
        var I = (d * h**3) / 12.0
        var K_l = mat.E * A / L
        var K_h = mat.E * I / L
        
        var step_size = (theta_end - theta_start) / Float64(num_steps - 1)
        
        for i in range(num_steps):
            var t_deg = theta_start + Float64(i) * step_size
            thetas_deg.append(t_deg)
            var t_rad_sweep = t_deg * pi / 180.0
            
            var factor: Float64 = 1.0
            if use_single_sided:
                factor = 2.0
            
            if i % 20 == 0:
                print("Step", i, "/", num_steps)
            
            var geom_sweep = get_rhombic_geometry(L, h, L0, h0, d, t_deg, L_flex)
            
            # (a) Ratio
            var r_val = calculate_amplification_ratio(mat, geom_sweep, use_single_sided)
            R_model.append(r_val)
            
            # (b) Input Stiffness
            var kin_val = calculate_input_stiffness(mat, geom_sweep, use_single_sided)
            K_in_model.append(kin_val / 1.0e6) # Convert N/m -> N/um
            
            # Frequencies
            var freqs_sweep = calculate_natural_frequencies(mat, geom_sweep, m_out, J_out, m_in)
            if len(freqs_sweep) >= 1: f1_model.append(freqs_sweep[0])
            else: f1_model.append(0.0)
            
            if len(freqs_sweep) >= 2: f2_model.append(freqs_sweep[1])
            else: f2_model.append(0.0)
            
            if len(freqs_sweep) >= 3: f3_model.append(freqs_sweep[2])
            else: f3_model.append(0.0)
            
            # Comparison Models 
            # 2. Eq 19 (Ref 23)
            var num19 = K_l * L**2 * sin(t_rad_sweep) * cos(t_rad_sweep)
            var den19 = 2.0 * K_h * cos(t_rad_sweep)**2 + K_l * L**2 * sin(t_rad_sweep)**2
            R_23.append((num19 / den19) * factor)
            
            # 3. Eq 20 (Ref 28)
            var cot_t = 1.0 / tan(t_rad_sweep)
            var tan_t = tan(t_rad_sweep)
            var term1 = (2.0 * K_h * cot_t) / (K_l * L)
            var term2 = L * tan_t
            R_28.append((L / (term1 + term2)) * factor)
            
            # 4. Eq 21 (Ref 31)
            var den21 = 12.0 * K_h * cos(t_rad_sweep)**2 + K_l * L**2 * sin(t_rad_sweep)**2
            R_31.append((num19 / den21) * factor)

        var fig = plt.figure()
        fig.set_size_inches(8, 12)
        
        # Subplot (a): Ratio
        var ax1 = plt.subplot(3, 1, 1)
        ax1.plot(thetas_deg, R_model, "r-")
        ax1.plot(thetas_deg, R_23, "k-")
        ax1.plot(thetas_deg, R_28, "g--")
        ax1.plot(thetas_deg, R_31, "b-.")
        ax1.set_ylabel("Amplification Ratio")
        
        var legend_a = builtins.list()
        legend_a.append("Proposed Model")
        legend_a.append("Ref [23]")
        legend_a.append("Ref [28]")
        legend_a.append("Ref [31]")
        ax1.legend(legend_a)
        
        ax1.grid(True)
        ax1.set_title("(a) Displacement amplification ratio")
        
        # Subplot (b): Input Stiffness
        var ax2 = plt.subplot(3, 1, 2)
        ax2.plot(thetas_deg, K_in_model, "r-")
        ax2.set_ylabel("Input Stiffness (N/um)")
        ax2.grid(True)
        ax2.set_title("(b) Input stiffness")
        
        # Subplot (c): Natural Frequencies
        var ax3 = plt.subplot(3, 1, 3)
        ax3.plot(thetas_deg, f1_model, "r-")
        ax3.plot(thetas_deg, f2_model, "g--")
        ax3.plot(thetas_deg, f3_model, "b-.")
        ax3.set_ylabel("Natural Frequency (Hz)")
        ax3.set_xlabel("Angle (deg)")
        
        var legend_c = builtins.list()
        legend_c.append("1st Mode")
        legend_c.append("2nd Mode")
        legend_c.append("3rd Mode")
        ax3.legend(legend_c)
        
        ax3.grid(True)
        ax3.set_title("(c) The first three in-plane natural frequencies")
        
        plt.tight_layout()
        plt.savefig("fig4_rhombic.png")
        print("Plot saved to fig4_rhombic.png")

    except e:
        print("Error with Python plot:", e)

fn get_parallel_geometry(
    L1: Float64, L2: Float64,
    h1: Float64, h2: Float64,
    L0: Float64, h0: Float64,
    d: Float64,
    H: Float64
) -> GeometricParameters:
    # See explanation in previous step.
    
    # Beam 1 (Input Port Half)
    var L_beam1 = sqrt((L0/2)**2 + (h0/2)**2)
    var h_beam1 = h0 
    var theta_beam1 = atan2(L0/2, h0/2)

    # Beam 2 (Flexure 1 - Horizontal)
    var L_beam2 = L1
    var h_beam2 = h1
    var theta_beam2 = 0.0

    # Beam 3 (Rigid Link/Flexure 2 - Angled)
    var L_beam3 = sqrt(L2**2 + H**2)
    var h_beam3 = h2
    var theta_beam3 = atan2(H, L2)

    # Beam 4 (Flexure Limb 4? - Horizontal)
    var L_beam4 = L1
    var h_beam4 = h1
    var theta_beam4 = 0.0

    return GeometricParameters(
        L_1=L_beam1, L_2=L_beam2, L_3=L_beam3, L_4=L_beam4,
        h_1=h_beam1, h_2=h_beam2, h_3=h_beam3, h_4=h_beam4,
        theta_1=theta_beam1, theta_2=theta_beam2, theta_3=theta_beam3, theta_4=theta_beam4,
        H=0.0, d=d
    )

fn get_aligned_geometry(
    L1: Float64, L2: Float64,
    h1: Float64, h2: Float64,
    L0: Float64, h0: Float64,
    d: Float64,
    theta_deg: Float64
) -> GeometricParameters:
    # Aligned-Type (Section 3.3)
    # Similar to Rhombic but with distinct flexure (h1) and link (h2) thicknesses/lengths.
    # All limb elements are aligned at angle theta.
    
    var theta_rad = theta_deg * pi / 180.0
    
    # Beam 1 (Input Port Half)
    var L_beam1 = sqrt((L0/2)**2 + (h0/2)**2)
    var h_beam1 = h0 
    var theta_beam1 = atan2(L0/2, h0/2)

    # Beam 2 (Flexure 1)
    var L_beam2 = L1
    var h_beam2 = h1
    var theta_beam2 = theta_rad

    # Beam 3 (Rigid Link / Flexure 2)
    var L_beam3 = L2
    var h_beam3 = h2
    var theta_beam3 = theta_rad

    # Beam 4 (Flexure 3 ? Check structure)
    # Parallel/Aligned structure is usually Flexure-Link-Flexure.
    # So L_beam2 = L1, L_beam3 = L2, L_beam4 = L1.
    var L_beam4 = L1
    var h_beam4 = h1
    var theta_beam4 = theta_rad

    return GeometricParameters(
        L_1=L_beam1, L_2=L_beam2, L_3=L_beam3, L_4=L_beam4,
        h_1=h_beam1, h_2=h_beam2, h_3=h_beam3, h_4=h_beam4,
        theta_1=theta_beam1, theta_2=theta_beam2, theta_3=theta_beam3, theta_4=theta_beam4,
        H=0.0, d=d
    )

fn solve_aligned_amplifier(
    mat: MaterialProperties,
    L1: Float64, L2: Float64,
    h1: Float64, h2: Float64,
    L0: Float64, h0: Float64,
    d: Float64,
    theta_deg: Float64,
    m_out: Float64,
    J_out: Float64,
    use_single_sided: Bool
) -> AmplifierResults:
    var geom = get_aligned_geometry(L1, L2, h1, h2, L0, h0, d, theta_deg)
    
    var res = AmplifierResults()
    res.amplification_ratio = calculate_amplification_ratio(mat, geom, use_single_sided)
    res.input_stiffness = calculate_input_stiffness(mat, geom, use_single_sided)
    res.natural_frequencies = calculate_natural_frequencies(mat, geom, m_out, J_out, 0.0)
    
    return res^

fn run_aligned_sweep(
    mat: MaterialProperties,
    L1: Float64, L2: Float64,
    h1: Float64, h2: Float64,
    L0_in: Float64, h0_in: Float64,
    d: Float64,
    m_out: Float64, J_out: Float64,
    use_single_sided: Bool,
    num_steps: Int,
    theta_start: Float64,
    theta_end: Float64
):
    print("\n------------------------------------------------")
    print("Generating Figure 6 (Aligned Amplifier) Sweep...")
    
    try:
        var plt = Python.import_module("matplotlib.pyplot")
        var builtins = Python.import_module("builtins")
        
        var thetas_deg = builtins.list()
        var R_model = builtins.list()
        var K_in_model = builtins.list()
        var f1_model = builtins.list()
        var f2_model = builtins.list()
        var f3_model = builtins.list()
        
        var step_size = (theta_end - theta_start) / Float64(num_steps - 1)
        
        for i in range(num_steps):
            var t_deg = theta_start + Float64(i) * step_size
            thetas_deg.append(t_deg)
            
            var geom = get_aligned_geometry(L1, L2, h1, h2, L0_in, h0_in, d, t_deg)
            
            # (a) Ratio
            var r_val = calculate_amplification_ratio(mat, geom, use_single_sided)
            R_model.append(r_val)
            
            # (b) Input Stiffness
            var kin_val = calculate_input_stiffness(mat, geom, use_single_sided)
            K_in_model.append(kin_val / 1.0e6) # N/um
            
            # (c) Frequencies
            var freqs_sweep = calculate_natural_frequencies(mat, geom, m_out, J_out, 0.0)
            if len(freqs_sweep) >= 1: f1_model.append(freqs_sweep[0])
            else: f1_model.append(0.0)
            
            if len(freqs_sweep) >= 2: f2_model.append(freqs_sweep[1])
            else: f2_model.append(0.0)
            
            if len(freqs_sweep) >= 3: f3_model.append(freqs_sweep[2])
            else: f3_model.append(0.0)
            
        var fig = plt.figure()
        fig.set_size_inches(8, 12)
        
        # Subplot (a): Ratio
        var ax1 = plt.subplot(3, 1, 1)
        ax1.plot(thetas_deg, R_model, "r-")
        ax1.set_ylabel("Amplification Ratio")
        ax1.set_title("(a) Displacement amplification ratio")
        ax1.grid(True)
        var legend_a = builtins.list()
        legend_a.append("Model")
        ax1.legend(legend_a)
        
        # Subplot (b): Input Stiffness
        var ax2 = plt.subplot(3, 1, 2)
        ax2.plot(thetas_deg, K_in_model, "r-")
        ax2.set_ylabel("Input Stiffness (N/um)")
        ax2.set_title("(b) Input stiffness")
        ax2.grid(True)
        
        # Subplot (c): Natural Frequencies
        var ax3 = plt.subplot(3, 1, 3)
        ax3.plot(thetas_deg, f1_model, "r-")
        ax3.plot(thetas_deg, f2_model, "g--")
        ax3.plot(thetas_deg, f3_model, "b-.")
        ax3.set_ylabel("Natural Frequency (Hz)")
        ax3.set_xlabel("Angle (deg)")
        ax3.set_title("(c) Natural frequencies")
        
        var legend_c = builtins.list()
        legend_c.append("1st Mode")
        legend_c.append("2nd Mode")
        legend_c.append("3rd Mode")
        ax3.legend(legend_c)
        ax3.grid(True)
        
        plt.tight_layout()
        plt.savefig("fig6_aligned.png")
        print("Saved fig6_aligned.png")

    except e:
        print("Error:", e)

fn solve_parallel_amplifier(
    mat: MaterialProperties,
    L1: Float64, L2: Float64,
    h1: Float64, h2: Float64,
    L0: Float64, h0: Float64,
    d: Float64,
    H: Float64,
    m_out: Float64,
    J_out: Float64,
    use_single_sided: Bool
) -> AmplifierResults:
    var geom = get_parallel_geometry(L1, L2, h1, h2, L0, h0, d, H)
    
    var res = AmplifierResults()
    res.amplification_ratio = calculate_amplification_ratio(mat, geom, use_single_sided)
    res.input_stiffness = calculate_input_stiffness(mat, geom, use_single_sided)
    res.natural_frequencies = calculate_natural_frequencies(mat, geom, m_out, J_out, 0.0)
    
    return res^

fn run_parallel_sweep(
    mat: MaterialProperties,
    L1: Float64, L2: Float64,
    h1: Float64, h2: Float64,
    L0_in: Float64, h0_in: Float64,
    d: Float64,
    m_out: Float64, J_out: Float64,
    use_single_sided: Bool,
    num_steps: Int,
    H_start: Float64,
    H_end: Float64
):
    print("\n------------------------------------------------")
    print("Generating Figure 5 (Parallel Amplifier) Sweep...")
    
    try:
        var plt = Python.import_module("matplotlib.pyplot")
        var builtins = Python.import_module("builtins")
        
        var H_vals = builtins.list()
        var R_model = builtins.list()
        var R_23 = builtins.list() # Eq 19
        var R_25 = builtins.list() # Eq 22
        var R_26 = builtins.list() # Eq 23
        var R_31 = builtins.list() # Eq 24
        
        var K_in_model = builtins.list()
        var f1_model = builtins.list()
        var f2_model = builtins.list()
        var f3_model = builtins.list()
        
        var step = (H_end - H_start) / Float64(num_steps - 1)
        
        for i in range(num_steps):
            var H = H_start + Float64(i) * step
            H_vals.append(H * 1000.0) # mm
            
            var geom = get_parallel_geometry(L1, L2, h1, h2, L0_in, h0_in, d, H)
            
            # (a) Ratio
            var r_val = calculate_amplification_ratio(mat, geom, use_single_sided)
            R_model.append(r_val)
            
            # (b) Input Stiffness
            var kin_val = calculate_input_stiffness(mat, geom, use_single_sided)
            K_in_model.append(kin_val / 1.0e6) # Convert N/m -> N/um
            
            # (c) Frequencies
            var freqs_sweep = calculate_natural_frequencies(mat, geom, m_out, J_out, 0.0) 
            if len(freqs_sweep) >= 1: f1_model.append(freqs_sweep[0])
            else: f1_model.append(0.0)
            
            if len(freqs_sweep) >= 2: f2_model.append(freqs_sweep[1])
            else: f2_model.append(0.0)
            
            if len(freqs_sweep) >= 3: f3_model.append(freqs_sweep[2])
            else: f3_model.append(0.0)
            
            # Comparsions
            var L0_total = L1 + L2
            # Paper Line 550: "h = atan(H/L0)" for reference models
            var theta_rad_ref = atan2(H, L0_total)
            var theta_deg_ref = theta_rad_ref * 180.0 / pi
            
            # Ref 23 (Eq 19)
            # Uses L (total length?) -> L0_total
            var val_23 = calculate_ratio_ref23(mat, L0_total, h1, d, theta_deg_ref)
            if use_single_sided: val_23 *= 2.0
            R_23.append(val_23)

            # Ref 25 (Eq 22)
            var val_25 = calculate_ratio_ref25_parallel(mat, L0_total, h1, d, H)
            if use_single_sided: val_25 *= 2.0
            R_25.append(val_25)
            
            # Ref 26 (Eq 23)
            var val_26 = calculate_ratio_ref26_parallel(mat, L0_total, h1, d, H)
            if use_single_sided: val_26 *= 2.0
            R_26.append(val_26)
            
            # Ref 31 (Eq 24)
            # Using best-effort reconstruction
            var val_31 = calculate_ratio_ref31_parallel(mat, L0_total, h1, h2, d, H, L2)
            if use_single_sided: val_31 *= 2.0
            R_31.append(val_31)
            
        var fig = plt.figure()
        fig.set_size_inches(8, 12)
        
        # Subplot (a): Ratio
        var ax1 = plt.subplot(3, 1, 1)
        ax1.plot(H_vals, R_model, "r-")
        ax1.plot(H_vals, R_23, "k-")
        ax1.plot(H_vals, R_25, "c--")
        ax1.plot(H_vals, R_26, "b-.")
        ax1.plot(H_vals, R_31, "m:")
        
        ax1.set_ylabel("Disp Ratio")
        ax1.set_title("(a) Displacement amplification ratio")
        
        var legend_a = builtins.list()
        legend_a.append("Model")
        legend_a.append("Ref 23")
        legend_a.append("Ref 25")
        legend_a.append("Ref 26")
        legend_a.append("Ref 31")
        ax1.legend(legend_a)
        ax1.grid(True)
        
        # Subplot (b): Input Stiffness
        var ax2 = plt.subplot(3, 1, 2)
        ax2.plot(H_vals, K_in_model, "r-")
        ax2.set_ylabel("Input Stiffness (N/um)")
        ax2.set_title("(b) Input stiffness")
        ax2.grid(True)
        
        # Subplot (c): Natural Frequencies
        var ax3 = plt.subplot(3, 1, 3)
        ax3.plot(H_vals, f1_model, "r-")
        ax3.plot(H_vals, f2_model, "g--")
        ax3.plot(H_vals, f3_model, "b-.")
        ax3.set_ylabel("Natural Frequency (Hz)")
        ax3.set_xlabel("H (mm)")
        ax3.set_title("(c) Natural frequencies")
        
        var legend_c = builtins.list()
        legend_c.append("1st Mode")
        legend_c.append("2nd Mode")
        legend_c.append("3rd Mode")
        ax3.legend(legend_c)
        ax3.grid(True)
        
        plt.tight_layout()
        plt.savefig("fig5_parallel.png")
        print("Saved fig5_parallel.png")

    except e:
        print("Error:", e)


fn main():
    var mat = MaterialProperties(E=71e9, rho=2770.0)
    var use_single_sided = True

    # ----------------------------------------------------------------
    # Rhombic Amplifier Config (Paper Section 3.1)
    # ----------------------------------------------------------------
    var L = 15.0e-3
    var h = 1.0e-3
    var L0 = 10.0e-3
    var h0 = 5.0e-3
    var d = 10.0e-3
    var theta_deg = 10.0
    var L_flex = 3.0e-3
    
    var volume_out = 4.0e-3 * 6.0e-3 * 10.0e-3
    var m_out = volume_out * mat.rho
    var J_out = (m_out / 12.0) * ((4.0e-3)**2 + (6.0e-3)**2)
    var m_in = 0.0

    print("Running Modular Rhombic Amplifier Analysis...")
    print("Parameters: theta =", theta_deg, "deg, h =", h*1e3, "mm, L =", L*1e3, "mm")

    var results = solve_rhombic_amplifier(
        mat, L, h, L0, h0, d, theta_deg, L_flex,
        m_out, J_out, m_in, use_single_sided
    )
    
    results.dump()

    var num_steps = 500
    var theta_start = 0.01
    var theta_end = 30.0
    run_rhombic_sweep(mat, L, h, L0, h0, d, L_flex, m_out, J_out, m_in, use_single_sided, num_steps, theta_start, theta_end)

    # ----------------------------------------------------------------
    # Parallel Amplifier Config (Paper Section 3.2)
    # ----------------------------------------------------------------
    print("\n------------------------------------------------")
    print("Running Modular Parallel Amplifier Analysis...")
    
    var L1_par = 3.0e-3
    var L2_par = 9.0e-3
    var h1_par = 0.5e-3
    var h2_par = 4.0e-3
    var L0_in_par = 10.0e-3
    var h0_in_par = 5.0e-3
    var d_par = 10.0e-3

    # For parallel, we explicitly calculate mass/inertia (Line 520)
    # m = 4x6x10 mm block. 
    var vol_par = 4.0e-3 * 6.0e-3 * 10.0e-3
    var m_out_par = vol_par * mat.rho
    var J_out_par = (m_out_par / 12.0) * ((4.0e-3)**2 + (6.0e-3)**2)
    # m_in neglected as per paper usually, or same as m_out? Paper says "m = ..." (singular lumped mass).
    # Assuming m corresponds to Output Mass.

    # 1. Single Point Analysis at theta = 10 deg (implied H)
    var theta_par_target = 10.0
    var theta_par_rad = theta_par_target * pi / 180.0
    # theta = atan(H/L2) => tan(theta) = H/L2 => H = L2 * tan(theta)
    var H_target = L2_par * tan(theta_par_rad)
    
    print("Parameters: equivalent theta =", theta_par_target, "deg, H =", H_target*1e3, "mm")
    
    var results_par = solve_parallel_amplifier(
        mat, L1_par, L2_par, h1_par, h2_par, L0_in_par, h0_in_par, d_par, 
        H_target, m_out_par, J_out_par, use_single_sided
    )
    
    results_par.dump()
    
    # 2. Sweep
    var H_start = 0.01e-3
    var H_end = 5.0e-3
    var num_steps_par = 500
    
    run_parallel_sweep(mat, L1_par, L2_par, h1_par, h2_par, L0_in_par, h0_in_par, d_par, m_out_par, J_out_par, use_single_sided, num_steps_par, H_start, H_end)

    # ----------------------------------------------------------------
    # Aligned Amplifier Config (Paper Section 3.3)
    # ----------------------------------------------------------------
    print("\n------------------------------------------------")
    print("Running Modular Aligned Amplifier Analysis...")
    
    var L1_ali = 3.0e-3
    var L2_ali = 9.0e-3
    var h1_ali = 0.5e-3
    var h2_ali = 4.0e-3
    var L0_in_ali = 10.0e-3
    var h0_in_ali = 5.0e-3
    var d_ali = 10.0e-3
    
    # Mass same as Parallel
    var m_out_ali = m_out_par
    var J_out_ali = J_out_par # same block
    
    # Sweep theta (h)
    var theta_start_ali = 0.01
    var theta_end_ali = 30.0
    var num_steps_ali = 100
    
    run_aligned_sweep(mat, L1_ali, L2_ali, h1_ali, h2_ali, L0_in_ali, h0_in_ali, d_ali, m_out_ali, J_out_ali, use_single_sided, num_steps_ali, theta_start_ali, theta_end_ali)