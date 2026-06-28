# two_analytical_ik_compare.cpp 修改说明

> 文件路径: `7DOF/6_28/two_analytical_ik_compare.cpp`

---

## 1. 文件头部改动

- 添加 `#define _USE_MATH_DEFINES`（Windows MinGW 兼容 `M_PI`）
- 添加 `#include <cstring>`（`memcpy` 使用）

---

## 2. 新增函数（PY 移植版）

在 `analytical_ik_paper_with_arm_angle_cal` 之后、`clampArmAngle` 之前插入约 240 行，完整移植了 `tsinghua_paper.py` 的几何求解链路。所有函数以 `_py` 结尾。

### solve_theta4_from_triangle_py
余弦定理求 θ₄: |SW|² = l_se² + l_ew² - 2·l_se·l_ew·cosθ₄

### compute_swe_from_target_py
从目标位姿求肩心 S、腕心 W、θ₄ 绝对值和 SW 单位向量

### elbow_from_arm_angle_py
臂角 θ₀ → 肘位置 E（圆几何：圆心 C + 半径 r × (cosθ₀·e₁ + sinθ₀·e₂)）

### solve_q123_from_swe_py
E、W、θ₄ 反解 q₁q₂q₃（遍历 s₂ 正负得两组 shoulder 解）

### extract_567_from_T47_py
从 T₄₇ 提取 wrist q₅q₆q₇（atan2 直接读矩阵元，s₆ 遍历±得两组）

### ik_one_arm_angle_py
给定单臂角的完整 IK：S→W→θ₄→E→q₁₂₃→FK到joint4→T₄₇→q₅₆₇

### analytical_ik_py
暴力扫描 [-π, π) 步长 0.01 rad（~628 次 IK 调用）选最优解

---

## 3. main 函数改动

原 `main` 只调 `select_optimal_ik_golden` 跑 C++ 版并打印结果。

新 `main` 改为**同场对比**：

1. 跑 `select_optimal_ik_golden`（C++ 原始版，golden section 搜索）
2. 跑 `analytical_ik_py`（PY 移植版，暴力扫描）
3. 对两版结果分别：
   - FK 验证 → 输出位置/姿态误差
   - 输出关节角
   - 记录耗时
4. 最终输出对比总结行：
   - C++ 版耗时 / PY 版耗时
   - 速度比 = PY 耗时 / C++ 耗时
   - 最大关节角偏差（度）

---

## 4. 两版关键差异

| 维度 | C++ 原始版 | PY 移植版 |
|------|-----------|----------|
| 臂角搜索 | golden section 15 次 | 暴力 ~628 次扫描 |
| wrist 偏置 | `sqrt(d6² + d7²)` 合并偏移 | 只减 `d7`，遗漏 `d6 = -2.5mm` |
| q₄ 符号 | `-acos(...)`（DH offset） | `+acos(...)`（py 原逻辑） |
| shoulder 解 | 从 R₀₃ ZYZ 分解 | 从 E 位置代入 FK 方程 |
| 求解耗时 | ~μs 级 | ~ms 级 |

> ⚠ PY 移植版严格按 `tsinghua_paper.py` 的数学逻辑翻译，**有意保留了 py 与 cpp 之间的原始差异**（如手腕偏置遗漏、q₄ 符号），以便验证和测量两套公式的实际差异。

---

## 5. 文件位置

- 修改后文件: `7DOF/6_28/two_analytical_ik_compare.cpp`
- 备份: `7DOF/6_28/two_analytical_ik_compare.cpp.bak`

---

## 6. WSL 编译运行

```bash
cd /mnt/c/Users/Lenovo/Desktop/7dof/7DOF/6_28

# 确保已安装 Eigen3
sudo apt install libeigen3-dev

# 编译
g++ -std=c++17 -O2 -o two_analytical_ik_compare \
    two_analytical_ik_compare.cpp ../r7_fk_lib.cpp -lEigen3

# 运行
./two_analytical_ik_compare
```
