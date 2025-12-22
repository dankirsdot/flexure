# Two-Port Dynamic Stiffness Model for Flexure Displacement Amplifier

This project will implement the general two-port dynamic stiffness model for three bridge-type flexure displacement amplifiers (Rhombic, Parallel, Aligned) as described in *Ling (2019)*.

I selected Mojo for education purposes, to get some idea about the language and its ecosystem. For now this file is an AI-generated plan with some manual edits for myself, with some basic structure to guide the implementation.

## Modules

### 1. `model.mojo`
This file will contain the physics and mathematics.
- **Structs**: To hold geometric parameters ($h, L, \theta, etc.$) and material properties ($E, \rho$).
- **Functions**:
    - `get_beam_stiffness(params, freq)`: Calculates $D_e(\omega)$ (Eq 1, 2).
    - `transform_matrix(local_matrix, angle)`: Rotates matrices (Eq 3, 5).
    - `get_transfer_matrix(stiffness_matrix)`: Converts stiffness to transfer matrix (Eq 11).
    - `condense_limb( ... )`: Combines multiple beam elements into one equivalent stiffness matrix (Eq 13).
    - `assemble_system( ... )`: Constructs the final system matrix (Eq 15).
    - `solve_...()`: Functions to find determinants (frequencies) and inverse (static/dynamic response).

### 2. `main.mojo`
This file drives the simulation.
- Define specific dimensions for the Rhombic/Parallel/Aligned amplifiers.
- Call `model` functions to build the system.
- Iterate through frequencies to find roots (natural frequencies).
- Solve for steady-state response at $\omega=0$ (static amplification/stiffness).

## Implementation Plan

### Step 1: Data Structures
**Goal**: Create a way to pass parameters around easily.
- Define a `struct` for `GeometricParameters` (Beam dimensions $L, h, b$, Angle $\theta$).
- Define a `struct` or `class` for `MaterialProperties` (Young's modulus `E`, density `rho`).

### Step 2: Beam Dynamic Stiffness Matrix
**Goal**: Implement Equation (1) and (2).
- Create a function `calc_beam_stiffness(geom, mat, omega)` returning a 6x6 Matrix.
- **Key Point**: The matrix elements ($\delta_1 .. \delta_8$) depend on frequency dependent parameters $\alpha$ and $\beta$.
- *Tip*: Mojo's SIMD or Tensor types can be useful for matrices, or use a simple 2D array/pointer wrapper if external linear algebra libs are not yet available.

### Step 3: Coordinate Transformation
**Goal**: Implement Equation (3) and (5).
- Implement the Rotation Matrix $R_i$ (Eq 5).
- Create a function to apply $R^T D_e R$.

### Step 4: Transfer Matrix & Condensation
**Goal**: Implement simplified limb logic (Eq 11, 12, 13).
- Implement conversion from Stiffness $D$ to Transfer $T$ (Eq 11).
- Implement matrix multiplication for $T_{total} = T_4 T_3 T_2 T_1$.
- Convert back from $T_{total}$ to equivalent Stiffness $D_{eq}$ (Eq 13).

### Step 5: System Assembly
**Goal**: Build the final matrix (Eq 15).
- Combine the condensed limb matrices into the global system matrix.
- Note the structure of Equation 15 – it's a sparse block matrix assembly.

### Step 6: Solvers
**Goal**: Extract results.
- **Statics ($\omega=0$)**:
    - Solve $F = K X$. Invert the matrix (or use Gaussian elimination).
    - Calculate Amplification Ratio $R$, Input Stiffness $K_{in}$ (Eq 18).
- **Dynamics ($\omega > 0$)**:
    - **Natural Frequency**: Find $\omega$ where $\det(D(\omega)) = 0$. You can use a root-finding algorithm (e.g., Bisection or Newton-Raphson) on the determinant function.

## Verification
We will verify the implementation using the **Rhombic Amplifier** case from Section 3.1 of the paper.
- **Parameters**: $h=1mm, L=15mm, \theta=10^\circ$, etc.
- **Expected Result**: Compare against Figure 4 curves.
