# 臂角可行区间计算 — 调试记录

## 背景

文件 `get_arm_angle_distrct.cpp` 实现了 7-DOF 冗余机械臂的**臂角（psi）可行区间**计算。
给定目标位姿 `T_target`，通过各关节的运动学约束，反解出 psi 的可行集合。

---

## 一、各关节类型

机械臂使用 MDH 参数，R₀₃ 按 ZY'Z'' 欧拉角分解，R₄₇ 同结构分解：

| 关节 | 提取公式 | 类型 |
|------|---------|------|
| J0 (θ₁) | `atan2(-R₀₃(1,2), -R₀₃(0,2))` | **atan2 型** |
| J1 (θ₂) | `acos(R₀₃(2,2))` | **acos 型** |
| J2 (θ₃) | `atan2(-R₀₃(2,1), R₀₃(2,0))` | **atan2 型** |
| J3 (θ₄) | 由目标位置唯一确定 | **常数**（不参与 psi 约束） |
| J4 (θ₅) | `atan2(-R₄₇(1,2), -R₄₇(0,2))` | **atan2 型** |
| J5 (θ₆) | `acos(R₄₇(2,2))` | **acos 型** |
| J6 (θ₇) | `atan2(-R₄₇(2,1), R₄₇(2,0))` | **atan2 型** |

---

## 二、区间分类逻辑

### atan2 型：显式 Δ 判断

`theta = atan2(num(ψ), den(ψ))`，num/den 均为 psi 的 sinusoid：

```
Delta = at² + bt² - ct²
  at = Bd·Cn - Bn·Cd
  bt = An·Cd - Ad·Cn
  ct = An·Bd - Ad·Bn
```

| Delta | 类型 | 处理 |
|-------|------|------|
| ≥ 0 | 循环型 | 返回全范围 `[-π, π]` |
| < 0 | 单调型 | 反解 psi 边界 |

### acos 型：三种情况显式分支

`theta = acos(f(ψ))`，`f(ψ) = A·cos + B·sin + C`，值域 `[C-r, C+r]`：

```
情况1 不可行：  f值域与限位完全不相交 → 返回空集
情况2 循环型：  f值域完全在限位范围内  → 返回全范围（对应 atan2 的 Delta≥0）
情况3 单调型：  部分重叠               → 反解 psi 边界（两段对称弧）
```

---

## 三、边界验证说明

每个可行区间有两个端点，但**并非全部来自关节限位**：

| 端点来源 | 含义 | 验证结果 |
|---------|------|---------|
| 关节到达限位的 psi | 真实约束边界 | 某关节恰好 = ±175° |
| ±π（wrap point） | `[-π,π)` 表示法的人工切断 | 无关节在限位 |

判断是否为 wrap point：
```cpp
bool is_wrap = (fabs(fabs(psi) - M_PI) < 1e-9);
```

---

## 四、发现并修复的 Bug

### Bug 1：`solve_psi` 使用 `tan` 导致 π 周期歧义

**原因**：`tan` 以 π 为周期，`tan(-175°) == tan(5°)`，solver 可能找到错误象限的解。

**原代码**：
```cpp
double psi_L = solve_psi(tan(theta_min));  // ← tan 有 π 周期歧义
// 验证时也用 tan 比较
double t1 = tan(th), t2 = tan_t;
return std::abs(t1-t2) < 1e-6;            // ← 同样有歧义
```

**修复**：改用 `cos/sin` 构建线性方程，`atan2` 验证：
```cpp
// theta = atan2(num, den) 等价条件: num*cos(t) - den*sin(t) = 0
double A = An*cos(theta_t) - Ad*sin(theta_t);
double B = Bn*cos(theta_t) - Bd*sin(theta_t);
double C = Cn*cos(theta_t) - Cd*sin(theta_t);

// 验证：直接比较 atan2 结果
double th = atan2(num, den);
double diff = std::abs(th - theta_t);
return diff < 1e-6 || std::abs(diff - 2*M_PI) < 1e-6;
```

### Bug 2：theta1 行列索引错误

**原因**：iv1 使用了 `(1,1)/(0,1)` 列，但 `eval_at_psi` 中 theta1 实际从 `(1,2)/(0,2)` 列提取。

**修复**：
```cpp
// 修复前（错误）
auto iv1 = atan2_joint_interval(
    -A_s(1,1), -B_s(1,1), -C_s(1,1),
    -A_s(0,1), -B_s(0,1), -C_s(0,1), ...);

// 修复后（正确）
auto iv1 = atan2_joint_interval(
    -A_s(1,2), -B_s(1,2), -C_s(1,2),
    -A_s(0,2), -B_s(0,2), -C_s(0,2), ...);
```

---

## 五、最终验证结果

目标位姿下各关节限位区间与最终可行区间：

```
===== 各关节ψ可行区间 =====
  theta1  [atan2] 全范围
  theta2  [acos ] 全范围
  theta3  [atan2] [4.16°, 180°] ∪ [-180°, -4.16°]   ← J2 限位
  theta5  [atan2] [85.02°, 180°] ∪ [-180°, 77.93°]   ← J4 限位
  theta6  [acos ] 全范围
  theta7  [atan2] 全范围

===== 最终可行区间 =====
  [85.02°, 180.00°]  — 上界 wrap，下界 J4=+175°
  [4.16°,  77.93°]   — 上界 J4=-175°，下界 J2=-175°
  [-180.00°, -4.16°] — 上界 J2=+175°，下界 wrap
```

边界验证（有效限位边界）：

| psi | 触发关节 |
|-----|---------|
| +85.02° | J4 = +175.0° ✅ |
| +4.16° | J2 = -175.0° ✅ |
| +77.93° | J4 = -175.0° ✅ |
| -4.16° | J2 = +175.0° ✅ |
| ±180° | wrap point（非限位） |
