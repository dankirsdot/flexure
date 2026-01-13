from math import pi, sin, cos, sinh, cosh, sqrt, tan, atan2

@fieldwise_init
struct GeometricParameters(Copyable):
    var L_1: Float64
    var L_2: Float64
    var L_3: Float64
    var L_4: Float64
    var h_1: Float64
    var h_2: Float64
    var h_3: Float64
    var h_4: Float64
    var theta_1: Float64
    var theta_2: Float64
    var theta_3: Float64
    var theta_4: Float64
    var H: Float64
    var d: Float64

    fn dump(self):
        print("L_1: ", self.L_1)
        print("L_2: ", self.L_2)
        print("L_3: ", self.L_3)
        print("L_4: ", self.L_4)
        print("h_1: ", self.h_1)
        print("h_2: ", self.h_2)
        print("h_3: ", self.h_3)
        print("h_4: ", self.h_4)
        print("theta_1: ", self.theta_1)
        print("theta_2: ", self.theta_2)
        print("theta_3: ", self.theta_3)
        print("theta_4: ", self.theta_4)
        print("H: ", self.H)
        print("d: ", self.d)

    fn area(self, h: Float64) -> Float64:
        return h * self.d

    fn moment_of_inertia(self, h: Float64) -> Float64:
        return (self.d * h**3) / 12.0

    fn clone(self) -> GeometricParameters:
        return GeometricParameters(
            L_1=self.L_1,
            L_2=self.L_2,
            L_3=self.L_3,
            L_4=self.L_4,
            h_1=self.h_1,
            h_2=self.h_2,
            h_3=self.h_3,
            h_4=self.h_4,
            theta_1=self.theta_1,
            theta_2=self.theta_2,
            theta_3=self.theta_3,
            theta_4=self.theta_4,
            H=self.H,
            d=self.d,
        )

@fieldwise_init
struct MaterialProperties(Copyable):
    var E: Float64
    var rho: Float64

    fn dump(self):
        print("E: ", self.E)
        print("rho: ", self.rho)

struct Matrix[size: Int](Copyable):
    var data: List[Float64]

    fn __init__(out self):
        self.data = List[Float64](capacity=Self.size * Self.size)
        for _ in range(Self.size * Self.size):
            self.data.append(0.0)

    fn __init__(out self, var data: List[Float64]):
        self.data = data^

    fn __copyinit__(out self, existing: Self):
        self.data = List[Float64](capacity=Self.size * Self.size)
        for i in range(len(existing.data)):
            self.data.append(existing.data[i])

    fn __getitem__(self, row: Int, col: Int) -> Float64:
        return self.data[row * Self.size + col]

    fn __setitem__(mut self, row: Int, col: Int, val: Float64):
        self.data[row * Self.size + col] = val

    fn dump(self):
        for i in range(Self.size):
            var line = String("")
            for j in range(Self.size):
                line += String(self[i, j]) + "  "
            print(line)

    fn zero(mut self):
        for i in range(Self.size * Self.size):
            self.data[i] = 0.0

    fn clone(self) -> Matrix[Self.size]:
        var res = Matrix[Self.size]()
        for i in range(Self.size * Self.size):
            res.data[i] = self.data[i]
        return res^

    fn transpose(self) -> Matrix[Self.size]:
        var res = Matrix[Self.size]()
        for i in range(Self.size):
            for j in range(Self.size):
                res[j, i] = self[i, j]
        return res^

    fn __matmul__(self, other: Matrix[Self.size]) -> Matrix[Self.size]:
        var res = Matrix[Self.size]()
        for i in range(Self.size):
            for j in range(Self.size):
                var sum: Float64 = 0.0
                for k in range(Self.size):
                    sum += self[i, k] * other[k, j]
                res[i, j] = sum
        return res^

    fn __neg__(self) -> Matrix[Self.size]:
        var res = Matrix[Self.size]()
        for i in range(Self.size):
            for j in range(Self.size):
                res[i, j] = -self[i, j]
        return res^

    fn __add__(self, other: Matrix[Self.size]) -> Matrix[Self.size]:
        var res = Matrix[Self.size]()
        for i in range(Self.size):
            for j in range(Self.size):
                res[i, j] = self[i, j] + other[i, j]
        return res^

    fn __sub__(self, other: Matrix[Self.size]) -> Matrix[Self.size]:
        var res = Matrix[Self.size]()
        for i in range(Self.size):
            for j in range(Self.size):
                res[i, j] = self[i, j] - other[i, j]
        return res^

fn cot(x: Float64) -> Float64:
    return cos(x) / sin(x)

fn csc(x: Float64) -> Float64:
    return 1.0 / sin(x)


fn get_beam_stiffness(
    mat: MaterialProperties,
    geom: GeometricParameters,
    omega: Float64,
    l: Float64,
    h: Float64,
) -> Matrix[6]:
    var A = geom.area(h)
    var I = geom.moment_of_inertia(h)

    var d = List[Float64](capacity=8)
    for _ in range(8):
        d.append(0.0)

    if omega < 1e-6:
        # Static limits
        d[0] = (mat.E * A) / l
        d[4] = -(mat.E * A) / l
        d[1] = (12.0 * mat.E * I) / l**3
        d[2] = -(6.0 * mat.E * I) / l**2
        d[3] = (4.0 * mat.E * I) / l
        d[5] = -(12.0 * mat.E * I) / l**3
        d[6] = (6.0 * mat.E * I) / l**2
        d[7] = (2.0 * mat.E * I) / l
    else:
        var alpha_sq = (omega**2 * l**2 * mat.rho) / mat.E
        var alpha = sqrt(alpha_sq)

        var beta_4 = (omega**2 * l**4 * mat.rho * A) / (mat.E * I)
        var beta = sqrt(sqrt(beta_4))

        var R = 1.0 - cos(beta) * cosh(beta)

        if alpha < 1e-4:
            # Axial static limit for small alpha to avoid division by zero
            d[0] = (mat.E * A) / l
            d[4] = -(mat.E * A) / l
        else:
            d[0] = (mat.E * A * alpha * cot(alpha)) / l
            d[4] = -(mat.E * A * alpha * csc(alpha)) / l

        d[1] = (mat.E * I * beta**3 * (cos(beta) * sinh(beta) + sin(beta) * cosh(beta))) / (R * l**3)
        d[2] = -(mat.E * I * beta**2 * (sin(beta) * sinh(beta))) / (R * l**2)
        d[3] = (mat.E * I * beta * (sin(beta) * cosh(beta) - cos(beta) * sinh(beta))) / (R * l)
        d[5] = -(mat.E * I * beta**3 * (sin(beta) + sinh(beta))) / (R * l**3)
        d[6] = (mat.E * I * beta**2 * (cosh(beta) - cos(beta))) / (R * l**2)
        d[7] = (mat.E * I * beta * (sinh(beta) - sin(beta))) / (R * l)

    var K = Matrix[6]()
    K[0, 0] = d[0];  K[0, 3] = d[4]
    K[1, 1] = d[1];  K[1, 2] = -d[2]; K[1, 4] = d[5];  K[1, 5] = d[6]
    K[2, 1] = -d[2]; K[2, 2] = d[3];  K[2, 4] = -d[6];  K[2, 5] = d[7]
    K[3, 0] = d[4];  K[3, 3] = d[0]
    K[4, 1] = d[5];  K[4, 2] = -d[6];  K[4, 4] = d[1];  K[4, 5] = d[2]
    K[5, 1] = d[6]; K[5, 2] = d[7];  K[5, 4] = d[2];  K[5, 5] = d[3]

    return K^

fn get_rotation_matrix_z(theta: Float64) -> Matrix[3]:
    var R = Matrix[3]()
    var c = cos(theta)
    var s = sin(theta)

    R[0, 0] = c;  R[0, 1] = s;  R[0, 2] = 0.0
    R[1, 0] = -s; R[1, 1] = c;  R[1, 2] = 0.0
    R[2, 0] = 0.0; R[2, 1] = 0.0; R[2, 2] = 1.0

    return R^

fn get_rotation_matrix(theta: Float64) -> Matrix[6]:
    var R = Matrix[6]()
    var R_z = get_rotation_matrix_z(theta)

    set_submatrix_3x3(R, R_z, 0, 0)
    set_submatrix_3x3(R, R_z, 3, 3)

    return R^

fn transform_matrix(K_hat: Matrix[6], theta: Float64) -> Matrix[6]:
    var R = get_rotation_matrix(theta)
    return R.transpose() @ K_hat @ R^

fn inverse_2x2(M: Matrix[2]) -> Matrix[2]:
    var res = Matrix[2]()
    var det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    if abs(det) < 1e-50:
        print("Warning: Singular 2x2 matrix")
        return res^

    var invDet = 1.0 / det
    res[0, 0] = M[1, 1] * invDet
    res[0, 1] = -M[0, 1] * invDet
    res[1, 0] = -M[1, 0] * invDet
    res[1, 1] = M[0, 0] * invDet
    return res^

fn inverse_3x3(M: Matrix[3]) -> Matrix[3]:
    var res = Matrix[3]()
    var t1 = M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1])
    var t2 = M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0])
    var t3 = M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0])
    var det = t1 - t2 + t3

    if abs(det) < 1e-50:
        print("Warning: Singular 3x3 matrix")
        return res^

    var invDet = 1.0 / det

    res[0, 0] = (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]) * invDet
    res[0, 1] = (M[0, 2] * M[2, 1] - M[0, 1] * M[2, 2]) * invDet
    res[0, 2] = (M[0, 1] * M[1, 2] - M[0, 2] * M[1, 1]) * invDet

    res[1, 0] = (M[1, 2] * M[2, 0] - M[1, 0] * M[2, 2]) * invDet
    res[1, 1] = (M[0, 0] * M[2, 2] - M[0, 2] * M[2, 0]) * invDet
    res[1, 2] = (M[1, 0] * M[0, 2] - M[0, 0] * M[1, 2]) * invDet

    res[2, 0] = (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0]) * invDet
    res[2, 1] = (M[2, 0] * M[0, 1] - M[0, 0] * M[2, 1]) * invDet
    res[2, 2] = (M[0, 0] * M[1, 1] - M[1, 0] * M[0, 1]) * invDet

    return res^

fn get_submatrix_3x3(K: Matrix[6], row_start: Int, col_start: Int) -> Matrix[3]:
    var res = Matrix[3]()
    for i in range(3):
        for j in range(3):
            res[i, j] = K[row_start + i, col_start + j]
    return res^

fn set_submatrix_3x3(mut K: Matrix[6], block: Matrix[3], row_start: Int, col_start: Int):
    for i in range(3):
        for j in range(3):
            K[row_start + i, col_start + j] = block[i, j]

fn get_transfer_from_stiffness(K: Matrix[6]) -> Matrix[6]:
    var k11 = get_submatrix_3x3(K, 0, 0)
    var k12 = get_submatrix_3x3(K, 0, 3)
    var k21 = get_submatrix_3x3(K, 3, 0)
    var k22 = get_submatrix_3x3(K, 3, 3)

    var k12_inv = inverse_3x3(k12)

    # According to the paper's convention T maps [x_j; F_j] -> [x_k; -F_k].
    # The negative sign on F_k is required for direct chain multiplication:
    # it satisfies equilibrium F_j_next + F_k = 0 => F_j_next = -F_k,
    # ensuring the input vector of the next element [x; F_j_next] matches [x; -F_k].
    var t11 = -(k12_inv @ k11)
    var t21 = -(k21 + k22 @ t11)
    var t22 = -(k22 @ k12_inv)
    var t12 = k12_inv^

    var T = Matrix[6]()
    set_submatrix_3x3(T, t11, 0, 0)
    set_submatrix_3x3(T, t12, 0, 3)
    set_submatrix_3x3(T, t21, 3, 0)
    set_submatrix_3x3(T, t22, 3, 3)

    return T^

fn get_stiffness_from_transfer(T: Matrix[6]) -> Matrix[6]:
    var t11 = get_submatrix_3x3(T, 0, 0)
    var t12 = get_submatrix_3x3(T, 0, 3)
    var t21 = get_submatrix_3x3(T, 3, 0)
    var t22 = get_submatrix_3x3(T, 3, 3)

    var t12_inv = inverse_3x3(t12)

    # Convert back from T ([x; -F] convention) to K (standard F=Kx).
    # Since T outputs -F_k and K requires +F_k, we must negate the rows
    # corresponding to forces (2nd block row) when deriving K from T.
    # Note: Eq. 13 in the paper has a typo and misses this negative sign.
    var k11 = -(t12_inv @ t11)
    var k22 = -(t22 @ t12_inv)
    var k12 = t12_inv^
    var k21 = -(t21 + t22 @ k11)

    var K = Matrix[6]()
    set_submatrix_3x3(K, k11, 0, 0)
    set_submatrix_3x3(K, k12, 0, 3)
    set_submatrix_3x3(K, k21, 3, 0)
    set_submatrix_3x3(K, k22, 3, 3)

    return K^

fn assemble_system_matrix_limb(mat: MaterialProperties, geom: GeometricParameters, omega: Float64) -> Matrix[6]:
    # Beam 1
    var K1_hat = get_beam_stiffness(mat, geom, omega, geom.L_1, geom.h_1)
    var K1 = transform_matrix(K1_hat, geom.theta_1)
    var T1 = get_transfer_from_stiffness(K1)

    # Beam 2
    var K2_hat = get_beam_stiffness(mat, geom, omega, geom.L_2, geom.h_2)
    var K2 = transform_matrix(K2_hat, geom.theta_2)
    var T2 = get_transfer_from_stiffness(K2)

    # Beam 3
    var K3_hat = get_beam_stiffness(mat, geom, omega, geom.L_3, geom.h_3)
    var K3 = transform_matrix(K3_hat, geom.theta_3)
    var T3 = get_transfer_from_stiffness(K3)

    # Beam 4
    var K4_hat = get_beam_stiffness(mat, geom, omega, geom.L_4, geom.h_4)
    var K4 = transform_matrix(K4_hat, geom.theta_4)
    var T4 = get_transfer_from_stiffness(K4)

    var T_total = T4 @ T3 @ T2 @ T1

    var K_limb = get_stiffness_from_transfer(T_total)

    return K_limb^

fn set_submatrix_9x9(mut K: Matrix[9], block: Matrix[3], row_start: Int, col_start: Int):
    for i in range(3):
        for j in range(3):
            K[row_start + i, col_start + j] = block[i, j]


fn assemble_full_system_matrix(mat: MaterialProperties, geom: GeometricParameters, omega: Float64) -> Matrix[9]:
    # Limb 1 (top left, 1st quadrant)
    var geom1 = geom.clone()
    var K_limb1 = assemble_system_matrix_limb(mat, geom1, omega)

    # Limb 2 (top right, 2nd quadrant)
    var geom2 = geom.clone()
    geom2.theta_1 = pi - geom2.theta_1
    geom2.theta_2 = pi - geom2.theta_2
    geom2.theta_3 = pi - geom2.theta_3
    geom2.theta_4 = pi - geom2.theta_4
    var K_limb2 = assemble_system_matrix_limb(mat, geom2, omega)

    # Limbs 3 and 4 are defined in the opposite direction compared to the paper.
    # It allows using the same geom object for all 4 limbs.

    # Limb 3 (bottom right, 3rd quadrant)
    var geom3 = geom.clone()
    geom3.theta_1 = pi + geom3.theta_1
    geom3.theta_2 = pi + geom3.theta_2
    geom3.theta_3 = pi + geom3.theta_3
    geom3.theta_4 = pi + geom3.theta_4
    var K_limb3 = assemble_system_matrix_limb(mat, geom3, omega)

    # Limb 4 (bottom left, 4th quadrant)
    var geom4 = geom.clone()
    geom4.theta_1 = -geom4.theta_1
    geom4.theta_2 = -geom4.theta_2
    geom4.theta_3 = -geom4.theta_3
    geom4.theta_4 = -geom4.theta_4
    var K_limb4 = assemble_system_matrix_limb(mat, geom4, omega)

    var K_sys = Matrix[9]()

    var k1_11 = get_submatrix_3x3(K_limb1, 0, 0)
    var k1_12 = get_submatrix_3x3(K_limb1, 0, 3)
    var k1_21 = get_submatrix_3x3(K_limb1, 3, 0)
    var k1_22 = get_submatrix_3x3(K_limb1, 3, 3)

    var k2_11 = get_submatrix_3x3(K_limb2, 0, 0)
    var k2_12 = get_submatrix_3x3(K_limb2, 0, 3)
    var k2_21 = get_submatrix_3x3(K_limb2, 3, 0)
    var k2_22 = get_submatrix_3x3(K_limb2, 3, 3)

    var k3_11 = get_submatrix_3x3(K_limb3, 0, 0)
    var k4_11 = get_submatrix_3x3(K_limb4, 0, 0)

    # Node 0 (InR): L1_Start + L4_End
    set_submatrix_9x9(K_sys, k1_11 + k4_11, 0, 0)

    # Node 1 (InL): L2_Start + L3_End
    set_submatrix_9x9(K_sys, k2_11 + k3_11, 3, 3)

    # Node 2 (Out): L1_End + L2_End
    set_submatrix_9x9(K_sys, k1_22 + k2_22, 6, 6)

    # Couplings
    # InR -> Out (L1)
    set_submatrix_9x9(K_sys, k1_12, 0, 6)
    set_submatrix_9x9(K_sys, k1_21, 6, 0)

    # InL -> Out (L2)
    set_submatrix_9x9(K_sys, k2_12, 3, 6)
    set_submatrix_9x9(K_sys, k2_21, 6, 3)

    return K_sys^

fn solve_static_response_9x9(K_sys: Matrix[9], F: Matrix[9]) -> Matrix[9]:
    var A = K_sys.clone()
    var B = F.clone()
    var n = 9
    for i in range(n):
        var max_row = i
        var max_val = abs(A[i, i])
        for k in range(i + 1, n):
            if abs(A[k, i]) > max_val:
                max_val = abs(A[k, i])
                max_row = k
        for col in range(i, n):
            var temp = A[i, col]
            A[i, col] = A[max_row, col]
            A[max_row, col] = temp
        var temp_b = B[i, 0]
        B[i, 0] = B[max_row, 0]
        B[max_row, 0] = temp_b

        if abs(A[i, i]) < 1e-20:
            print("Singular 9x9")
            return Matrix[9]()

        for k in range(i + 1, n):
            var factor = A[k, i] / A[i, i]
            for j in range(i, n):
                A[k, j] = A[k, j] - factor * A[i, j]
            B[k, 0] = B[k, 0] - factor * B[i, 0]

    var X = Matrix[9]()
    for i in range(n - 1, -1, -1):
        var sum = 0.0
        for j in range(i + 1, n):
            sum += A[i, j] * X[j, 0]
        X[i, 0] = (B[i, 0] - sum) / A[i, i]
    return X^

fn calculate_amplification_ratio(mat: MaterialProperties, geom: GeometricParameters, use_single_sided: Bool = False) -> Float64:
    var K_sys = assemble_full_system_matrix(mat, geom, 0.0)

    var F = Matrix[9]()
    var f = 100.0
    F[0, 0] = f
    F[3, 0] = -f

    var X = solve_static_response_9x9(K_sys, F)

    var x_out_y = X[7, 0]
    var x_in1_x = X[0, 0]
    var x_in2_x = X[3, 0]

    if abs(x_in1_x - x_in2_x) < 1e-12:
        return 0.0

    # Differential Convention: Output / Total Input Stroke
    var ratio = x_out_y / (x_in1_x - x_in2_x)
    
    if use_single_sided:
        return ratio * 2.0
    
    return ratio

fn calculate_input_stiffness(mat: MaterialProperties, geom: GeometricParameters, use_single_sided: Bool = False) -> Float64:
    var K_sys = assemble_full_system_matrix(mat, geom, 0.0)

    var F = Matrix[9]()
    var f = 100.0
    F[0, 0] = f
    F[3, 0] = -f

    var X = solve_static_response_9x9(K_sys, F)

    var x_in1_x = X[0, 0]
    var x_in2_x = X[3, 0]

    if abs(x_in1_x - x_in2_x) < 1e-12:
        return 1e20 # Rigid

    # Differential Convention: Force / Total Input Stroke
    var k_in = f / (x_in1_x - x_in2_x)
    
    if use_single_sided:
        return k_in * 2.0
        
    return k_in

fn determinant_9x9(K: Matrix[9]) -> Float64:
    var A = K.clone()
    var n = 9
    var det = 1.0

    for i in range(n):
        var max_row = i
        var max_val = abs(A[i, i])
        for k in range(i + 1, n):
            if abs(A[k, i]) > max_val:
                max_val = abs(A[k, i])
                max_row = k

        if max_row != i:
            for col in range(i, n):
                var temp = A[i, col]
                A[i, col] = A[max_row, col]
                A[max_row, col] = temp
            det = -det

        if abs(A[i, i]) < 1e-50:
            return 0.0

        det *= A[i, i]

        for k in range(i + 1, n):
            var factor = A[k, i] / A[i, i]
            for j in range(i, n):
                A[k, j] = A[k, j] - factor * A[i, j]

    return det

fn calculate_natural_frequencies(
    mat: MaterialProperties,
    geom: GeometricParameters,
    m_out: Float64,
    J_out: Float64,
    m_in: Float64,
) -> List[Float64]:
    var start_f = 100.0
    var end_f = 10000.0
    var step_f = 10.0

    var frequencies = List[Float64]()
    var prev_det = 0.0
    var first_run = True

    var f = start_f
    while f < end_f:
        if len(frequencies) >= 3:
            break

        var omega = f * 2.0 * pi

        var K_sys = assemble_full_system_matrix(mat, geom, omega)

        var inertia_out = -(omega**2) * m_out
        K_sys[6, 6] = K_sys[6, 6] + inertia_out
        K_sys[7, 7] = K_sys[7, 7] + inertia_out

        var rot_inertia_out = -(omega**2) * J_out
        K_sys[8, 8] = K_sys[8, 8] + rot_inertia_out

        var inertia_in = -(omega**2) * m_in
        K_sys[0, 0] += inertia_in
        K_sys[1, 1] += inertia_in
        K_sys[3, 3] += inertia_in
        K_sys[4, 4] += inertia_in

        var det = determinant_9x9(K_sys)

        if not first_run:
            if det * prev_det < 0.0:
                # Zero crossing detected
                var f_prev = f - step_f
                # Linear interpolation for root
                var f_root = f_prev - prev_det * step_f / (det - prev_det)
                frequencies.append(f_root)
            elif abs(det) < 1e-12:
                frequencies.append(f)

        prev_det = det
        first_run = False
        f += step_f

    return frequencies^

fn calculate_ratio_ref23(
    mat: MaterialProperties,
    L: Float64,
    h: Float64,
    d: Float64,
    theta_deg: Float64,
) -> Float64:
    var theta_rad = theta_deg * pi / 180.0
    var A = d * h
    var I = (d * h**3) / 12.0
    var K_l = mat.E * A / L
    var K_h = mat.E * I / L

    var num = K_l * L**2 * sin(theta_rad) * cos(theta_rad)
    var den = 2.0 * K_h * cos(theta_rad)**2 + K_l * L**2 * sin(theta_rad)**2
    return num / den

fn calculate_ratio_ref28(
    mat: MaterialProperties,
    L: Float64,
    h: Float64,
    d: Float64,
    theta_deg: Float64,
) -> Float64:
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


fn calculate_ratio_ref31(
    mat: MaterialProperties,
    L: Float64,
    h: Float64,
    d: Float64,
    theta_deg: Float64,
) -> Float64:
    var theta_rad = theta_deg * pi / 180.0
    var A = d * h
    var I = (d * h**3) / 12.0
    var K_l = mat.E * A / L
    var K_h = mat.E * I / L

    var num = K_l * L**2 * sin(theta_rad) * cos(theta_rad)
    var den = 12.0 * K_h * cos(theta_rad)**2 + K_l * L**2 * sin(theta_rad)**2
    return num / den

fn calculate_ratio_ref25_parallel(
    mat: MaterialProperties,
    L0: Float64,
    h_flex: Float64,
    d: Float64,
    H: Float64,
) -> Float64:
    # Eq 22
    # R = (K_l L_0^2 sin h cos^3 h) / (2K_h + K_l L_0^2 cos^2 h sin h)
    
    var theta_rad = atan2(H, L0) 
    
    var A = d * h_flex
    var I = (d * h_flex**3) / 12.0
    var K_l = mat.E * A / L0
    var K_h = mat.E * I / L0
    
    var s = sin(theta_rad)
    var c = cos(theta_rad)
    
    var num = K_l * L0**2 * s * c**3
    var den = 2.0 * K_h + K_l * L0**2 * c**2 * s
    
    return num / den

fn calculate_ratio_ref26_parallel(
    mat: MaterialProperties,
    L0: Float64,
    h_flex: Float64,
    d: Float64,
    H: Float64,
) -> Float64:
    # Eq 23
    # R = (K_l L0 H) / (4 K_h + K_l H^2)
    
    var A = d * h_flex
    var I = (d * h_flex**3) / 12.0
    var K_l = mat.E * A / L0
    var K_h = mat.E * I / L0
    
    var num = K_l * L0 * H
    var den = 4.0 * K_h + K_l * H**2
    
    return num / den

fn calculate_ratio_ref31_parallel(
    mat: MaterialProperties,
    L0: Float64,
    h_flex: Float64,
    h_link: Float64,
    d: Float64,
    H: Float64,
    L_link: Float64
) -> Float64:
    # Eq 24 (Ref 31)
    # Attempting to reconstruct from fragmented OCR text.
    
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