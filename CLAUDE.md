# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Run

All root-level C++ files use C++17 and depend on Eigen3. Build by linking the shared FK library as needed:

```bash
# Simple standalone files (no FK dependency)
g++ -std=c++17 -o quatarnion_test quatarnion_test.cpp
g++ -std=c++17 -o test test.cpp -lEigen3

# Files using forward_kinematics — must link r7_fk_lib.cpp
g++ -std=c++17 -o numerical_ik_r7 numerical_ik_r7.cpp r7_fk_lib.cpp
g++ -std=c++17 -o fk_ik_loop fk_ik_loop.cpp r7_fk_lib.cpp
g++ -std=c++17 -o check_lie_algebra_rot_error check_lie_algebra_rot_error.cpp r7_fk_lib.cpp -lEigen3

# The two large analytical IK + cycling test binaries
g++ -std=c++17 -O2 -o analytical_ik analytical_ik.cpp r7_fk_lib.cpp -lEigen3
g++ -std=c++17 -O2 -o ik_cycling_test ik_cycling_test.cpp r7_fk_lib.cpp -lEigen3

# Pure FK (r7_fk.cpp has its own FK, unlike the others that use r7_fk_lib.cpp)
g++ -std=c++17 -o r7_fk r7_fk.cpp

# Standalone helpers
g++ -std=c++17 -o transform_to_matrices transform_to_matrices.cpp
g++ -std=c++17 -o fk_ik_loop fk_ik_loop.cpp r7_fk_lib.cpp

# Python
python3 triangle_angle.py

# Sigma7Device library
cd Sigma7Device && /usr/local/Qt-5.15.2/bin/qmake Sigma7Device.pro && make
```

## Architecture

This repo contains two unrelated pieces:

### 1. 7-DOF robot arm kinematics (repo root)

Kinematics for a 7-DOF redundant serial manipulator. **All angles are in degrees, all positions in mm.** Uses Modified DH (MDH) convention.

**Shared FK library — `r7_fk_lib.cpp`:**
The canonical forward kinematics implementation. Provides `forward_kinematics(const double q_deg[7], double T[4][4], int joint=7)` where `joint` selects which link's pose (1-7, default 7 = end-effector). Most other files declare this as `extern` and link against it. The DH table is hardcoded here and duplicated in several files.

**DH parameters (MDH):**

| i | αᵢ₋₁ (°) | aᵢ₋₁ (mm) | dᵢ (mm) |
|---|-----------|------------|---------|
| 1 | 0 | 0 | 84+95 |
| 2 | -90 | 0 | 0 |
| 3 | 90 | 0 | 215+260 |
| 4 | -90 | 0 | 0 |
| 5 | 90 | 0 | 415+60 |
| 6 | -90 | 0 | -2.5 |
| 7 | 90 | 0 | 145 |

Zero-configuration end-effector position: (0, -2.5, 1274) mm.

**Key kinematics files:**

- **`analytical_ik.cpp`** (~1730 lines) — The main analytical IK implementation based on the paper "Analytical Inverse Kinematic Computation for 7-DOF Redundant Manipulators With Joint Limits and Its Application to Redundancy Resolution." Core function: `analytical_ik_paper(T_target, q_init, psi, q_out)` — solves IK for a given target pose, initial guess, and arm angle ψ. Also includes `select_optimal_ik()` (exhaustive ψ sampling) and `select_optimal_ik_golden()` (golden-section search over ψ). Joint limits: ±175° (≈±3.054 rad) on all joints. Arm angle feasible intervals: [0, 0.9021] and [1.1019, 2.6463] rad.

- **`ik_cycling_test.cpp`** (~1820 lines) — Near-duplicate of `analytical_ik.cpp` with additional cycling/validation logic. Tests IK consistency across configuration changes. Uses the same `analytical_ik_paper()`, `select_optimal_ik_golden()`, `score_solution()` functions.

- **`golden_ik.h`** — Shared golden-section search helpers for arm angle optimization. Depends on `analytical_ik_paper()` and `score_solution()` being available (does not define them).

- **`numerical_ik_r7.cpp`** — Numerical IK using damped least squares (Levenberg-Marquardt) with nullspace projection for redundancy resolution. Interface: `numerical_ik(T_target, q_init, q_out, max_iter, pos_tol, rot_tol, verbose, q_prefer)`. Uses finite-difference Jacobian (dq=1e-4 rad), position rows scaled by 1/100 to match rotation rows, line search over {1.0, 0.5, 0.25} step sizes. Converges in 3-7 iterations, ~28-48 μs.

- **`r7_fk.cpp`** — Standalone FK binary with its own FK implementation (not the shared lib). Interactive, accepts 7 joint angles and prints end-effector pose (position + RPY + quaternion).

- **`fk_ik_loop.cpp`** — FK→IK closed-loop validation: user inputs target joint angles, FK computes pose, IK solves back, compares.

- **`analytical_ik_test.cpp`** — Minimal test scaffold for the analytical IK (hardcoded target pose, prints θ₄).

- **`check_lie_algebra_rot_error.cpp`** — Compares two rotation error formulations: axis-angle from relative rotation matrix vs Lie algebra log-SO(3). Uses Eigen.

**Other root-level utilities:**
- **`quatarnion_test.cpp`** — Quaternion to Euler angles (roll/pitch/yaw, ZYX convention).
- **`transform_to_matrices.cpp`** — Extracts 8 column-major 4x4 homogeneous transforms from a `double[8][16]` array, reports rotation axis/angle and translation.
- **`test.cpp`** — Minimal arm angle clamping test.
- **`triangle_angle.py`** — Law of cosines: given three sides, compute triangle angles.

**Important conventions:**
- Joint angle units: **degrees** everywhere (FK input, IK input/output). Convert internally to radians where needed.
- Position units: **mm** in FK output and IK target poses.
- MDH convention: α in degrees in the DH table, converted to radians inside `mdh_transform()`.
- `analytical_ik.cpp` and `ik_cycling_test.cpp` share ~90% of their code. Changes to IK logic likely need to be applied to both.
- The FK library `r7_fk_lib.cpp` has `-DIK_LIB_ONLY` mode that excludes `main()` for use as a compilation unit in other binaries.

### 2. Sigma7Device (Sigma7Device/)

A Qt/C++ shared library (`libSigma7Device.so`) for controlling Force Dimension haptic devices (Sigma.7, Omega, Lambda) used in surgical robotics. Built with qmake (Qt 5.15), no Qt GUI modules (Qt is used for the build system only).

**Key dependencies:**
- Force Dimension DHD/DRD SDK (bundled in `dhd/` — headers in `dhd/include/`, libs for Linux x64 and Windows)
- Ruckig motion library (external, at `../3rdParty/lib/ruckig/`)
- libusb-1.0, pthread

**Core classes and their roles:**
- `SigmaDeviceImpl` (Sigma7DeviceImpl.h) — main device interface. Implements `IOmegaDevice` to manage up to 2 haptic arms. Handles device init, force feedback, position locking, wrist alignment, clutch detection, filtering (OneEuro/MovingAverage), and coordinate transforms.
- `LambdaNavigate` / `LambdaViscosity` / `LambdaDataCalcul` — lambda-specific kinematics, viscosity force computation, and encoder-to-pose calculations.
- `myOmega` / `myOmegaWrist` / `myOmegaWrist4` — Omega device forward kinematics (delta robot geometry) and wrist kinematics variants.
- `ConeProjectionCalculator` — joint limit enforcement via cone projection.
- `trajectorygenerator` — trajectory generation (via ruckig).
- `PositionAlgorithm` — position computation utilities.
- `lambda.h` — shared data definitions: structs for raw/calculated lambda data, wrist joint enums, scaling factors, arm ID constants.

**Key constants (from headers):**
- `LEFT_OMEGA=0`, `RIGHT_OMEGA=1` — arm indices
- `ARM_ID_0`–`ARM_ID_3` — arm-to-lambda assignments
- Wrist joints: `FIRST_JOINT=yaw`, `SECOND_JOINT=pitch`, `THIRD_JOINT=roll`
