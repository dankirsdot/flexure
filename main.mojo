from model import *
from python import Python
from math import pi, sin, cos, tan, sqrt, atan2

fn main():
    var mat = MaterialProperties(E=71e9, rho=2770.0)
    
    # parameters as per paper
    # var L_0_p = 10.0e-3
    # var h_0_p = 5.0e-3
    # var L_1_p = 10.0e-3
    # var h_1_p = 5.0e-3
    # var L_2_p = 9.0e-3
    # var h_2_p = 1.0e-3
    # var H_p = 1.0e-3
    # var d_p = 1.0e-3

    # parameters as per model
    # var L_1 = sqrt((L_0_p/2)**2 + (h_0_p/2)**2)
    # var h_1 = h_0_p # treat te block as very stiff short beam
    # var theta_1 = atan2(L_0_p, h_0_p)
    
    # var L_2 = L_1_p
    # var h_2 = h_1_p
    # var theta_2 = pi - 10.0 * pi / 180.0

    # var L_3 = sqrt(L_2_p**2 + H_p**2)
    # var h_3 = h_2_p
    # var theta_3 = pi - 10.0 * pi / 180.0

    # var L_4 = L_2
    # var h_4 = h_2
    # var theta_4 = theta_2

    # ----------------------------------------------------------------
    # Verification Case: Rhombic Amplifier
    # ----------------------------------------------------------------
    # Paper Parameters (Section 3.1)
    var L_p = 15.0e-3
    var h_p = 1.0e-3
    var L0_p = 10.0e-3
    var h0_p = 5.0e-3
    var d_p = 10.0e-3
    var H_p = 0.0

    var theta_deg = 10.0
    var theta_rad = theta_deg * pi / 180.0

    # Flexure Split (L_p = L2 + L3 + L4). Uniform beam.
    # L2 = 3mm, L3 = 9mm, L4 = 3mm.
    var L_flex_p = 3.0e-3
    var L_rigid_p = L_p - 2.0 * L_flex_p

    # Beam 1 (Input Port Half)
    var L_1 = sqrt((L0_p/2)**2 + (h0_p/2)**2)
    var h_1 = h0_p 
    var theta_1 = atan2(L0_p/2, h0_p/2)

    # Beam 2 (Flexure 1)
    var L_2 = L_flex_p
    var h_2 = h_p
    var theta_2 = theta_rad

    # Beam 3 (Rigid Link)
    var L_3 = L_rigid_p # sqrt(L_rigid_p**2 + H_p**2)
    var h_3 = h_p
    var theta_3 = theta_2 # atan2(H_p, L_rigid_p)

    # Beam 4 (Flexure 2)
    var L_4 = L_flex_p
    var h_4 = h_p
    var theta_4 = theta_2

    var params_rhombic = GeometricParameters(
        L_1=L_1, L_2=L_2, L_3=L_3, L_4=L_4,
        h_1=h_1, h_2=h_2, h_3=h_3, h_4=h_4,
        theta_1=theta_1, theta_2=theta_2, theta_3=theta_3, theta_4=theta_4,
        H=H_p, d=d_p
    )

    # Amplification Ratio
    print("\n[Statics] Calculating Amplification Ratio...")
    var R = calculate_amplification_ratio(mat, params_rhombic)
    print("Amplification Ratio R =", R)
    print("Expected (approx cot 10) =", 1.0/tan(theta_rad))
    
    # Natural Frequency
    var volume_out = 4.0e-3 * 6.0e-3 * 10.0e-3
    var mass_in = 0.0
    var mass_out = volume_out * mat.rho
    
    print("\n[Dynamics] Calculating First Natural Frequencies...")
    var freqs = calculate_natural_frequencies(mat, params_rhombic, mass_out, mass_in)
    print("Natural Frequencies:", freqs[0], freqs[1], freqs[2])
    
    # ------------------------------------------------
    # Figure 4 Replication
    # ------------------------------------------------
    print("\n------------------------------------------------")
    print("Generating Figure 4 (a, b, c)...")
    
    try:
        var plt = Python.import_module("matplotlib.pyplot")
        var np = Python.import_module("numpy")
        
        var builtins = Python.import_module("builtins")
        var thetas_deg = builtins.list()
        
        # Data for (a)
        var R_model = builtins.list()
        var R_23 = builtins.list()
        var R_28 = builtins.list()
        var R_31 = builtins.list()
        
        # Data for (b)
        var K_in_model = builtins.list()
        
        # Data for (c)
        var f1_model = builtins.list()
        var f2_model = builtins.list()
        var f3_model = builtins.list()
        
        # Theoretical Constants
        var L_p = 15.0e-3
        var h_p = 1.0e-3
        var A = d_p * h_p
        var I = (d_p * h_p**3) / 12.0
        var K_l = mat.E * A / L_p
        var K_h = mat.E * I / L_p
        
        var num_steps = 100
        var theta_start = 0.01
        var theta_end = 30.0
        var step_size = (theta_end - theta_start) / Float64(num_steps - 1)
        
        print("Starting sweep...")
        for i in range(num_steps):
            var t_deg = theta_start + Float64(i) * step_size
            thetas_deg.append(t_deg)
            var t_rad_sweep = t_deg * pi / 180.0
            
            if i % 10 == 0:
                print("Step", i, "/", num_steps)
            
            var geom_sweep = GeometricParameters(
                L_1=L_1, L_2=L_flex_p, L_3=L_rigid_p, L_4=L_flex_p,
                h_1=h_1, h_2=h_p, h_3=h_p, h_4=h_p,
                theta_1=theta_1, theta_2=t_rad_sweep, theta_3=t_rad_sweep, theta_4=t_rad_sweep,
                H=0.0, d=d_p
            )
            
            # (a) Ratio
            var r_val = calculate_amplification_ratio(mat, geom_sweep)
            R_model.append(r_val)
            
            # (b) Input Stiffness
            var kin_val = calculate_input_stiffness(mat, geom_sweep)
            K_in_model.append(kin_val / 1.0e6) # Convert N/m -> N/um
            
            # Frequencies
            var freqs_sweep = calculate_natural_frequencies(mat, geom_sweep, mass_out, mass_in)
            if len(freqs_sweep) >= 1: f1_model.append(freqs_sweep[0])
            else: f1_model.append(0.0)
            
            if len(freqs_sweep) >= 2: f2_model.append(freqs_sweep[1])
            else: f2_model.append(0.0)
            
            if len(freqs_sweep) >= 3: f3_model.append(freqs_sweep[2])
            else: f3_model.append(0.0)
            # 2. Eq 19 (Ref 23)
            # R = (Kl L^2 sin cos) / (2 Kh cos^2 + Kl L^2 sin^2)
            var num19 = K_l * L_p**2 * sin(t_rad_sweep) * cos(t_rad_sweep)
            var den19 = 2.0 * K_h * cos(t_rad_sweep)**2 + K_l * L_p**2 * sin(t_rad_sweep)**2
            R_23.append(num19 / den19)
            
            # 3. Eq 20 (Ref 28)
            # R = L / ( (2 Kh cot) / (Kl L) + L tan )
            var cot_t = 1.0 / tan(t_rad_sweep)
            var tan_t = tan(t_rad_sweep)
            var term1 = (2.0 * K_h * cot_t) / (K_l * L_p)
            var term2 = L_p * tan_t
            R_28.append(L_p / (term1 + term2))
            
            # 4. Eq 21 (Ref 31)
            # R = (Kl L^2 sin cos) / (12 Kh cos^2 + Kl L^2 sin^2)
            var den21 = 12.0 * K_h * cos(t_rad_sweep)**2 + K_l * L_p**2 * sin(t_rad_sweep)**2
            R_31.append(num19 / den21)
            
        var fig = plt.figure()
        fig.set_size_inches(8, 12)
        
        # Subplot (a): Ratio
        var ax1 = plt.subplot(3, 1, 1)
        ax1.plot(thetas_deg, R_model, "r-")
        ax1.plot(thetas_deg, R_23, "k-")
        ax1.plot(thetas_deg, R_28, "g--")
        ax1.plot(thetas_deg, R_31, "b-.")
        ax1.set_ylabel("Amplification Ratio")
        # ax1.set_xlabel("Angle (deg)")
        
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
        # ax2.set_xlabel("Angle (deg)")
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
        plt.savefig("fig4_replication.png")
        print("Plot saved to fig4_replication.png")
        
    except e:
        print("Error with Python plot:", e)