#include <iostream>
#include <vector>
#include <utility> 
#include <cmath>
#include <chrono>
#include <Eigen/Geometry>
#include <iomanip>

using std::vector;
using std::pair;
using namespace std;

#define CONTROL_HZ  500.0
#define DT          (1.0 / CONTROL_HZ)
#define N_JOINTS    7
#define N_PSI       7

//计算出的臂角可行区间
constexpr double ARM_ANGLE_DISTRICT_1_1 = 0.0;
constexpr double ARM_ANGLE_DISTRICT_1_2 = 0.9021;
constexpr double ARM_ANGLE_DISTRICT_2_1 = 1.1019;
constexpr double ARM_ANGLE_DISTRICT_2_2 = 2.646297;

// 关节限位（弧度）175*M_PI/180 = 3.054326
static const double Q_MIN[7] = {-3.054326, -3.054326, -3.054326, -3.054326, -3.054326, -3.054326, -3.054326};
static const double Q_MAX[7] = { 3.054326,  3.054326,  3.054326,  3.054326,  3.054326,  3.054326, 3.054326};

template <typename Derived>
void printEigenMatrix(const Eigen::MatrixBase<Derived>& mat,
                      const std::string& name = "matrix")
{
    std::cout << name << ":\n";
    std::cout << mat << std::endl;
}

const double dh[7][4] = {
    {   0,  0,   84+95,   0 },
    {  -90,  0,   0,       0 },
    {  90,  0,  215+260,   0 },
    {  -90,  0,   0,       0 },
    {  90,  0,  415+60,    0 },
    {  -90,  0, -2.5,       0 },
    {  90,  0,  145,        0 },
};
void forward_kinematics(const double joint_angles[7], double T_out[4][4], int joint = 7);
void fK_eigen(const double joint_angles[7], Eigen::Matrix3d& R, Eigen::Vector3d& p, int joint = 7);
void getPoseFromArray(const double T[4][4], Eigen::Matrix3d& R, Eigen::Vector3d& p)
{
    R << T[0][0], T[0][1], T[0][2],
         T[1][0], T[1][1], T[1][2],
         T[2][0], T[2][1], T[2][2];

    p << T[0][3], T[1][3], T[2][3];
}

// 角度差归一化到 [-180°, 180°]，处理 atan2 分支切割引起的 2π 跳变
static inline double unwrap_delta(double q, double q_ref) {
    double dq = q - q_ref;
    while (dq >  180.0) dq -= 360.0;
    while (dq < -180.0) dq += 360.0;
    return dq;
}



// 区间求交集工具
static vector<pair<double,double>> intersect(
    const vector<pair<double,double>>& a,
    const vector<pair<double,double>>& b)
{
    vector<pair<double,double>> res;
    for (auto& x : a)
        for (auto& y : b) {
            double lo = std::max(x.first, y.first);
            double hi = std::min(x.second, y.second);
            if (hi > lo + 1e-9) res.push_back({lo, hi});
        }
    return res;
}

// 全范围 [-pi, pi)
static vector<pair<double,double>> full_range() {
    return {{-M_PI, M_PI}};
}

// atan2型关节的可行区间
// theta(psi) = atan2(An*sin+Bn*cos+Cn, Ad*sin+Bd*cos+Cd)
// Delta = at^2 + bt^2 - ct^2
//   Delta>0: 循环型，theta可遍历全范围，不施加约束 → 返回full_range
//   Delta<0: 单调型，用theta_min/max反解psi边界
//   Delta=0: 奇异型，暂返回full_range（保守处理）
static vector<pair<double,double>> atan2_joint_interval(
    double An, double Bn, double Cn,
    double Ad, double Bd, double Cd,
    double theta_min, double theta_max)
{
    double at = Bd*Cn - Bn*Cd;
    double bt = An*Cd - Ad*Cn;
    double ct = An*Bd - Ad*Bn;
    double Delta = at*at + bt*bt - ct*ct;

    if (Delta >= 0.0) return full_range(); // 循环型或奇异型

    // 单调型：确定方向（在psi=0处导数符号）
    bool increasing = (bt + ct > 0.0);

    // 对给定 theta_target，解 psi：atan2(num,den)=theta_t
    // 条件: num*cos(t) - den*sin(t) = 0，即 A*sin+B*cos+C=0
    auto solve_psi = [&](double theta_t) -> double {
        double ct = cos(theta_t), st = sin(theta_t);
        double A = An*ct - Ad*st;
        double B = Bn*ct - Bd*st;
        double C = Cn*ct - Cd*st;
        double r = sqrt(A*A + B*B);
        if (std::abs(C) > r + 1e-9) return NAN;
        double phi = atan2(B, A);
        double sv = std::clamp(-C/r, -1.0, 1.0);
        double alpha = asin(sv);
        double p1 = alpha - phi;
        double p2 = M_PI - alpha - phi;
        auto wrap = [](double v) {
            while (v < -M_PI) v += 2*M_PI;
            while (v >= M_PI) v -= 2*M_PI;
            return v;
        };
        p1 = wrap(p1); p2 = wrap(p2);
        // 用 atan2 验证，避免 tan 的 π 周期歧义
        auto check = [&](double psi) {
            double num = An*sin(psi)+Bn*cos(psi)+Cn;
            double den = Ad*sin(psi)+Bd*cos(psi)+Cd;
            double th = atan2(num, den);
            double diff = std::abs(th - theta_t);
            return diff < 1e-6 || std::abs(diff - 2*M_PI) < 1e-6;
        };
        if (check(p1)) return p1;
        if (check(p2)) return p2;
        return p1;
    };

    double psi_L = solve_psi(theta_min);
    double psi_U = solve_psi(theta_max);
    if (std::isnan(psi_L) || std::isnan(psi_U)) return full_range();

    double psi_start = increasing ? psi_L : psi_U;
    double psi_end   = increasing ? psi_U : psi_L;

    if (psi_start <= psi_end)
        return {{psi_start, psi_end}};
    else
        return {{psi_start, M_PI}, {-M_PI, psi_end}};
}

// acos型关节：theta(psi) = acos(A*cos+B*sin+C)，theta in [0,pi]
// 极值点(theta最小/最大)对应的psi: 式(37)(38)
// 返回theta在[theta_min, theta_max]内可行的psi区间
static vector<pair<double,double>> acos_joint_interval(
    double A, double B, double C,
    double theta_min, double theta_max)
{
    // theta(psi) = acos(f(psi)), f(psi) = A*cos + B*sin + C
    // f 的值域: [f_min, f_max] = [C-r, C+r]
    double r = sqrt(A*A + B*B);
    double f_min = C - r;
    double f_max = C + r;

    // 关节限位对应的 f 边界（acos 单调递减：theta大→f小）
    double f_lo = cos(theta_max); // theta_max → f 下界
    double f_hi = cos(theta_min); // theta_min → f 上界

    // ---- 情况1：不可行（f 值域与限位要求完全不相交）----
    if (f_max < f_lo - 1e-9 || f_min > f_hi + 1e-9) return {};

    // ---- 情况2：循环型（f 值域完全在限位范围内，所有 psi 均可行）----
    // 对应 atan2 型的 Delta >= 0 情况
    if (f_min >= f_lo - 1e-9 && f_max <= f_hi + 1e-9) return full_range();

    // ---- 情况3：单调型（部分重叠，需反解 psi 边界）----
    double f_lo_eff = std::max(f_min, f_lo);
    double f_hi_eff = std::min(f_max, f_hi);

    double phi = atan2(B, A); // f(psi) 取最大值处的 psi
    // f >= f_lo_eff: |psi - phi| <= half_lo
    double half_lo = acos(std::clamp((f_lo_eff - C)/r, -1.0, 1.0));
    // f <= f_hi_eff: |psi - phi| >= half_hi
    double half_hi = acos(std::clamp((f_hi_eff - C)/r, -1.0, 1.0));

    auto wrap = [](double v) {
        while (v >= M_PI) v -= 2*M_PI;
        while (v < -M_PI) v += 2*M_PI;
        return v;
    };
    vector<pair<double,double>> res;
    auto add = [&](double lo, double hi) {
        lo = wrap(lo); hi = wrap(hi);
        if (lo <= hi) res.push_back({lo, hi});
        else { res.push_back({lo, M_PI}); res.push_back({-M_PI, hi}); }
    };
    // 可行区间: [phi-half_lo, phi-half_hi] ∪ [phi+half_hi, phi+half_lo]
    if (half_hi < half_lo - 1e-9) {
        add(phi - half_lo, phi - half_hi);
        add(phi + half_hi, phi + half_lo);
    } else {
        add(phi - half_lo, phi + half_lo);
    }
    return res;
}




// 计算两个齐次变换矩阵的差异
void compareTransforms(const double T1[4][4], const double T2[4][4]) {
    // 提取平移部分 (假设矩阵按行主序存储，T[row][col])
    double tx1 = T1[0][3], ty1 = T1[1][3], tz1 = T1[2][3];
    double tx2 = T2[0][3], ty2 = T2[1][3], tz2 = T2[2][3];
    
    // 位置误差：平移向量的欧几里得距离
    double dx = tx1 - tx2;
    double dy = ty1 - ty2;
    double dz = tz1 - tz2;
    double position_error = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // 提取旋转矩阵部分 (3x3)
    double R1[3][3], R2[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            R1[i][j] = T1[i][j];
            R2[i][j] = T2[i][j];
        }
    }
    
    // 计算相对旋转矩阵 ΔR = R1 * R2^T
    // 首先计算 R2 的转置 R2_T
    double R2_T[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            R2_T[i][j] = R2[j][i];
    
    // 矩阵乘法 ΔR = R1 * R2_T
    double delta_R[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                delta_R[i][j] += R1[i][k] * R2_T[k][j];
    
    // 通过迹+反对称部分计算旋转角度，用 atan2 对小角度更稳定
    double trace = delta_R[0][0] + delta_R[1][1] + delta_R[2][2];
    double s = std::sqrt((delta_R[2][1] - delta_R[1][2]) * (delta_R[2][1] - delta_R[1][2]) +
                         (delta_R[0][2] - delta_R[2][0]) * (delta_R[0][2] - delta_R[2][0]) +
                         (delta_R[1][0] - delta_R[0][1]) * (delta_R[1][0] - delta_R[0][1])) / 2.0;
    double c = (trace - 1.0) / 2.0;
    double angle_error_rad = std::atan2(s, c);  // 弧度
    double angle_error_deg = angle_error_rad * 180.0 / M_PI;  // 转换为度
    
    // 输出对比结果
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "位置误差 (Position error): " << position_error << " units\n";
    std::cout << "姿态误差 (Orientation error): " << angle_error_rad << " rad ("
              << angle_error_deg << " deg)\n";
}


//新建反对成矩阵
Eigen::Matrix3d skew(const Eigen::Vector3d& u)
{
    Eigen::Matrix3d ux;
    ux <<  0.0,   -u.z(),  u.y(),
           u.z(),  0.0,   -u.x(),
          -u.y(),  u.x(),  0.0;
    return ux;
}

void fK_eigen(const double joint_angles[7], Eigen::Matrix3d& R, Eigen::Vector3d& p, int joint) {
    double T[4][4];
    forward_kinematics(joint_angles, T, joint);
    getPoseFromArray(T, R, p);
}


/**
 * @brief 计算两个平面的夹角：
 *       参考平面：theta_3 = 0 时，shoulder(frame2)、elbow(frame3)、wrist(frame5) 确定的平面
 *       实际平面：theta_3 ≠ 0 时的同一平面
 * @param q 当前角度
 * @return psi 返回弧度
 */
double arm_plane_angle(const double q[7]) {
    // 1. 获取三个关键点Eigen::Matrix3d R_elbow; Eigen::Vector3d p_elbow;
    Eigen::Matrix3d R_elbow; 
    Eigen::Vector3d p_elbow;
    fK_eigen(q, R_elbow, p_elbow, 4);  // 肘关节位置

    double T[4][4] = {};
    forward_kinematics(q, T, 7);
    Eigen::Matrix3d R_ee; 
    Eigen::Vector3d p_ee;
    getPoseFromArray(T, R_ee, p_ee);

    // 2. 肩中心（固定）
    Eigen::Vector3d p_s(0, 0, dh[0][2]);

    // 3. 腕中心
    Eigen::Vector3d l_7_wt(0, 0, sqrt(pow(dh[6][2],2) + pow(dh[5][2],2)));
    Eigen::Vector3d p_w = p_ee - R_ee * l_7_wt;

    // 4. 肩→腕单位向量（旋转轴）
    Eigen::Vector3d sw = p_w - p_s;
    Eigen::Vector3d u = sw.normalized();

    // 5. 参考方向：世界Z轴投影到垂直u的平面
    Eigen::Vector3d z_world(0, 0, 1);
    Eigen::Vector3d v_ref = z_world - z_world.dot(u) * u;
    if (v_ref.norm() < 1e-6) {
        // sw近似平行Z轴时改用X轴
        Eigen::Vector3d x_world(1, 0, 0);
        v_ref = x_world - x_world.dot(u) * x_world;
    }
    v_ref.normalize();

    // 6. 当前肘向量投影到垂直u的平面
    Eigen::Vector3d se = p_elbow - p_s;
    Eigen::Vector3d v_e = se - se.dot(u) * u;
    if (v_e.norm() < 1e-6) return 0.0;  // 奇异：肘在肩腕连线上
    v_e.normalize();

    // 7. 带符号角度
    double cos_psi = std::clamp(v_ref.dot(v_e), -1.0, 1.0);
    double sin_psi = u.dot(v_ref.cross(v_e));
    return std::atan2(sin_psi, cos_psi);
}


/**
 * @brief 单调型 (Delta3 < 0) 求解臂角可行区间
 * 
 * 单调型特点：
 *   - θ(ψ) 在整个 [-π, π) 上单调递增或递减
 *   - 无极值点
 *   - 通过 θ_min/θ_max 直接解对应的 ψ，再根据单调方向确定区间
 * 
 * @param A_n, B_n, C_n    分子系数: num = A_n*sinψ + B_n*cosψ + C_n
 * @param A_d, B_d, C_d    分母系数: den = A_d*sinψ + B_d*cosψ + C_d
 * @param a_t, b_t, c_t    导数系数（用于确定单调方向）
 * @param theta_min        关节角下限（弧度）
 * @param theta_max        关节角上限（弧度）
 * @param intervals        输出：可行臂角区间列表
 * @return true  求解成功
 * @return false 求解失败（如无解）
 */
bool solveMonotonicIntervals(
    double A_n, double B_n, double C_n,
    double A_d, double B_d, double C_d,
    double a_t, double b_t, double c_t,
    double theta_min, double theta_max,
    vector<pair<double, double>>& intervals)
{
    intervals.clear();

    // ============================================================
    // 1. 确定单调方向
    //    在 ψ=0 处，导数的分子为: a_t*sin(0) + b_t*cos(0) + c_t = b_t + c_t
    //    分母恒正，所以导数的符号由 (b_t + c_t) 决定
    // ============================================================
    double sign_val = b_t + c_t;
    bool increasing = (sign_val > 0.0);   // true: θ随ψ递增; false: 递减

    // ============================================================
    // 2. 辅助函数：对于给定的 tan_target，解对应的 ψ（归一化到 [-π, π)）
    //    方程： (A_n sinψ + B_n cosψ + C_n) / (A_d sinψ + B_d cosψ + C_d) = tan_target
    // ============================================================
    auto solve_psi_for_tan = [&](double tan_target) -> double {
        double A = A_n - tan_target * A_d;
        double B = B_n - tan_target * B_d;
        double C = C_n - tan_target * C_d;

        double r = sqrt(A * A + B * B);
        if (std::abs(C) > r) return NAN;   // 无解

        // 解 A sinψ + B cosψ = -C
        double phi = atan2(B, A);           // A sinψ + B cosψ = r sin(ψ + phi)
        double sin_val = -C / r;
        double alpha = asin(sin_val);       // ψ + phi = alpha 或 π - alpha

        double psi_1 = alpha - phi;
        double psi_2 = M_PI - alpha - phi;

        // 归一化到 [-π, π)
        auto wrap = [](double x) {
            while (x < -M_PI) x += 2 * M_PI;
            while (x >=  M_PI) x -= 2 * M_PI;
            return x;
        };
        psi_1 = wrap(psi_1);
        psi_2 = wrap(psi_2);

        // 检验哪个 ψ 真的满足 tan(θ(ψ)) = tan_target
        auto check = [&](double psi) {
            double num = A_n * sin(psi) + B_n * cos(psi) + C_n;
            double den = A_d * sin(psi) + B_d * cos(psi) + C_d;
            double theta_here = atan2(num, den);
            return std::abs(tan(theta_here) - tan_target) < 1e-9;
        };

        if (check(psi_1)) return psi_1;
        if (check(psi_2)) return psi_2;

        // fallback
        return psi_1;
    };

    // ============================================================
    // 3. 对 theta_min 和 theta_max 求对应的臂角边界
    // ============================================================
    double psi_L = solve_psi_for_tan(tan(theta_min));
    double psi_U = solve_psi_for_tan(tan(theta_max));

    if (std::isnan(psi_L) || std::isnan(psi_U)) {
        // 无解：限位内没有可行臂角
        return false;
    }

    // ============================================================
    // 4. 根据单调方向确定可行区间
    // ============================================================
    double psi_start, psi_end;
    if (increasing) {
        psi_start = psi_L;
        psi_end   = psi_U;
    } else {
        psi_start = psi_U;
        psi_end   = psi_L;
    }

    // ============================================================
    // 5. 输出区间（考虑跨越 -π 的情况）
    // ============================================================
    if (psi_start <= psi_end) {
        intervals.push_back({psi_start, psi_end});
    } else {
        intervals.push_back({psi_start,  M_PI});
        intervals.push_back({-M_PI, psi_end});
    }

    return true;
}


enum sigualrity_type{
    safty = 0,
    shoulder_singular = -1,
    elbow_singular = -2,
    wrist_singular = -3,
    shoulder_wrist_singular = -4
};

/**
 * @brief 需要传入肩关节和腕关节之间的连线
 * @param q_init 为角度
 * @return sigualrity_type 奇异类型
 */
sigualrity_type check_near_singularity(double q_init[7], Eigen::Vector3d x_0_sw){
    //1、肩部奇异，条件为x_0_sw × Z0 = 0。 一般不会到达
    
    // 基座标系 Z 轴 (0,0,1)
    Eigen::Vector3d z0(0.0, 0.0, 1.0);
    // 计算叉积
    Eigen::Vector3d cross_result = x_0_sw.cross(z0);
    if(fabs(cross_result.norm()) < 0.05){ //当腕部距离z轴线 5cm
        return shoulder_singular;
    }
    //2、肘部奇异，计算出来的 theta4 不为 0/pi。 一般只有工作空间边界才有
    else if(fabs(q_init[3]) < 10.0 || fabs(q_init[3]+180) < 10.0 ||fabs(q_init[3]-180) < 10.0 ){
        return elbow_singular;
    }

    //3、腕部奇异 theta6 角度接近 0，关节 5、 6 共线的情况
    else if(fabs(q_init[5]) < 10.0 || fabs(q_init[5]+180) < 10.0 ||fabs(q_init[5]-180) < 10.0 ){
        return wrist_singular;
    }

    //4、文章中提及的 sin theta2，6 为零的情况
    else if(fabs(q_init[1]) < 10.0 || fabs(q_init[1]+180) < 10.0 ||fabs(q_init[1]-180) < 10.0 ){
        return shoulder_wrist_singular;
    }
    else{
        return safty;
    }

}


void analytical_ik_paper_with_arm_angle_cal(const double T_target[4][4]/* , const double q_init[7], double psi */, double q_out[7]){
    double arm_angle = 0.597493;

    Eigen::Matrix3d R_0_desire;
    Eigen::Vector3d P_0_desire;
    //获取base下的SW坐标,及其单位向量
    Eigen::Vector3d x_0_sw; // 在 base frame中 获取sw向量
    Eigen::Vector3d u_0_sw;//x_sw_0 的单位向量

    getPoseFromArray(T_target,R_0_desire,P_0_desire);//获取目标位姿的posi rot
    Eigen::Vector3d l_0_bs(0,0,dh[0][2]);
    Eigen::Vector3d l_7_wt(0,0,sqrt(pow(dh[6][2],2) + pow(dh[5][2],2)));
    x_0_sw = P_0_desire - l_0_bs - R_0_desire *l_7_wt;

    printf("x_0_sw的位置为: %f,%f,%f\n",x_0_sw(0),x_0_sw(1),x_0_sw(2));
    
    double norm_x_0_sw = x_0_sw.norm();
    if (norm_x_0_sw > 1e-12) {
        u_0_sw = x_0_sw / norm_x_0_sw;
    } else {
        u_0_sw.setZero();  
    }

    //==============解 theta 4=========== checked
    double cos_theta4 = (pow(norm_x_0_sw,2) - pow(dh[2][2],2) - pow(dh[4][2],2)) / ( 2* dh[2][2] * dh[4][2]);
    
    if (cos_theta4 < -1.0) cos_theta4 = -1.0;
    if (cos_theta4 > 1.0) cos_theta4 = 1.0;
    double theta_4 = -std::acos(cos_theta4);
    printf("------------第四个关节角度为：%f \n",theta_4*180/M_PI);
    q_out[3] = - theta_4 * 180.0 / M_PI;//负号为offset

    // ===== 构造 R_3_4 (theta_4 已求出，提到前面供后面复用) =====
    double c4 = std::cos(theta_4);
    double s4 = std::sin(theta_4);

    //暂时不使用fk计算（坐标系不一样），由于 j4是绕着 base frame 的 -y轴旋转的，按照旋转为：
    Eigen::Matrix3d R_3_4;
    R_3_4 <<  c4, 0.0, -s4,
            0.0, 1.0, 0.0,
            s4, 0.0,  c4;


    auto normalizeAngle = [](double angle) {
        while (angle > M_PI) angle -= 2.0 * M_PI;
        while (angle < -M_PI) angle += 2.0 * M_PI;
        return angle;
    };
    //第1个轴是绕base frame 的 z
    auto Rz = [](double theta) {
        double c = std::cos(theta);
        double s = std::sin(theta);
        Eigen::Matrix3d R;
        R << c, -s, 0.0,
             s,  c, 0.0,
             0.0, 0.0, 1.0;
        return R;
    };
    //第二个轴是绕base frame 的 -y
    auto Ry_neg = [](double theta) {
        double c = std::cos(theta);
        double s = std::sin(theta);
        Eigen::Matrix3d R;
        R <<  c, 0.0, -s,
             0.0, 1.0, 0.0,
             s, 0.0,  c;
        return R;
    };
    const double eps = 1e-9;

    //求当臂角为0时候的theta_1_ref plane 和theta_1_ref plane 对应公式 14
    double theta_1_ref = 0.0;
    double theta_2_ref = 0.0;

    // ================= 求 theta_1_ref, theta_2_ref =================
    //
    // 目标形式：
    // x_0_sw = Rz(theta_1_ref) * Ry(-theta_2_ref) * p
    //
    // p = R_2_3(theta3 = 0) * (l_3_se + R_3_4(theta4) * l_4_ew)
    Eigen::Matrix3d R_2_3_theta3_zero = Eigen::Matrix3d::Identity();
    Eigen::Vector3d l_3_se(0.0, 0.0, dh[2][2]);
    Eigen::Vector3d l_4_ew(0.0, 0.0, dh[4][2]);

    // 已知向量 p
    Eigen::Vector3d p = R_2_3_theta3_zero * (l_3_se + R_3_4 * l_4_ew);

    // 现在解：x_0_sw = Rz(theta_1_ref) * Ry(-theta_2_ref) * p
    // [ c1·c2,  -s1,  -c1·s2 ]                                   
    // [ s1·c2,   c1,  -s1·s2 ]                                   
    // [ s2,       0,     c2   ]   
    //有两个未知数，3组方程。直接解耦求解

    const double eps_1 = 1e-12;

    const double x = x_0_sw(0);
    const double y = x_0_sw(1);
    const double z = x_0_sw(2);

    const double px = p(0);
    const double py = p(1);
    const double pz = p(2);

    const double r = std::sqrt(px * px + pz * pz);

    if (r < eps_1) {
        printf("求 theta_2_ref 失败：p 在 y 轴附近，出现奇异。\n");
    } else {
        double cos_arg = z / r;

        if (cos_arg > 1.0) cos_arg = 1.0;
        if (cos_arg < -1.0) cos_arg = -1.0;

        const double alpha = std::atan2(px, pz);
        const double beta = std::acos(cos_arg);

        double theta_2_candidate_1 = alpha + beta;
        double theta_2_candidate_2 = alpha - beta;

        double theta_1_candidate_1 = 0.0;
        double theta_1_candidate_2 = 0.0;

        {
            double c2 = std::cos(theta_2_candidate_1);
            double s2 = std::sin(theta_2_candidate_1);

            double a = c2 * px - s2 * pz;

            theta_1_candidate_1 =
                std::atan2(y, x) - std::atan2(py, a);
        }

        {
            double c2 = std::cos(theta_2_candidate_2);
            double s2 = std::sin(theta_2_candidate_2);

            double a = c2 * px - s2 * pz;

            theta_1_candidate_2 =
                std::atan2(y, x) - std::atan2(py, a);
        }

        theta_1_candidate_1 = normalizeAngle(theta_1_candidate_1);
        theta_2_candidate_1 = normalizeAngle(theta_2_candidate_1);

        theta_1_candidate_2 = normalizeAngle(theta_1_candidate_2);
        theta_2_candidate_2 = normalizeAngle(theta_2_candidate_2);

        // 用正运动学残差选一个更匹配的解
        Eigen::Vector3d x_check_1 =
            Rz(theta_1_candidate_1) * Ry_neg(theta_2_candidate_1) * p;

        Eigen::Vector3d x_check_2 =
            Rz(theta_1_candidate_2) * Ry_neg(theta_2_candidate_2) * p;

        double err_1 = (x_check_1 - x_0_sw).norm();
        double err_2 = (x_check_2 - x_0_sw).norm();

        if (err_1 <= err_2) {
            theta_1_ref = theta_1_candidate_1;
            theta_2_ref = theta_2_candidate_1;
        } else {
            theta_1_ref = theta_1_candidate_2;
            theta_2_ref = theta_2_candidate_2;
        }

        printf("theta_1_ref candidate 1: %f deg, theta_2_ref candidate 1: %f deg, err: %.12f\n",
            theta_1_candidate_1 * 180.0 / M_PI,
            theta_2_candidate_1 * 180.0 / M_PI,
            err_1);

        printf("theta_1_ref candidate 2: %f deg, theta_2_ref candidate 2: %f deg, err: %.12f\n",
            theta_1_candidate_2 * 180.0 / M_PI,
            theta_2_candidate_2 * 180.0 / M_PI,
            err_2);

        printf("selected theta_1_ref: %f deg\n", theta_1_ref * 180.0 / M_PI);
        printf("selected theta_2_ref: %f deg\n", theta_2_ref * 180.0 / M_PI);
    }
    // ================= 根据 arm_angle 求 theta_1, theta_2, theta_3 =================
    //
    // 这里假设肩部前三个旋转为：
    // R_0_1 = Rz(theta_1)
    // R_1_2 = Ry(-theta_2)
    // R_2_3 = Rz(theta_3)
    //
    // 所以：
    // R_0_3 = Rz(theta_1) * Ry(-theta_2) * Rz(theta_3)

    double theta_1 = 0.0;
    double theta_2 = 0.0;
    double theta_3 = 0.0;

    // 1. 先构造 arm_angle = 0 时的参考 R_0_3_ref
    // theta_3_ref = 0
    Eigen::Matrix3d R_0_3_ref =
        Rz(theta_1_ref) * Ry_neg(theta_2_ref) * Rz(0.0);

    // 2. 构造 [u_0_sw x] 反对称矩阵
    Eigen::Matrix3d u_cross;
    u_cross << 0.0,        -u_0_sw(2),  u_0_sw(1),
            u_0_sw(2),   0.0,       -u_0_sw(0),
            -u_0_sw(1),   u_0_sw(0),  0.0;

    // 3. 根据论文公式构造 A_s, B_s, C_s
    Eigen::Matrix3d A_s = u_cross * R_0_3_ref;
    Eigen::Matrix3d B_s = -u_cross * u_cross * R_0_3_ref;
    Eigen::Matrix3d C_s = (u_0_sw * u_0_sw.transpose()) * R_0_3_ref;

    {
        //---------------------求 theta1的可行范围---------------------
        double A_n_1 = -A_s(1,1);
        double B_n_1 = -B_s(1,1);
        double C_n_1 = -C_s(1,1);
        double A_d_1 = -A_s(0,1);
        double B_d_1 = -B_s(0,1);
        double C_d_1 = -C_s(0,1);

        // 计算 a_t, b_t, c_t（公式27）
        double a_t1 = B_d_1 * C_n_1 - B_n_1 * C_d_1;
        double b_t1 = A_n_1 * C_d_1 - A_d_1 * C_n_1;
        double c_t1 = A_n_1 * B_d_1 - A_d_1 * B_n_1;

        // 判别式
        double Delta1 = a_t1*a_t1 + b_t1*b_t1 - c_t1*c_t1;

        //确认theta1为单调，按照关节范围对应最大最小
        double theta_min_1 = -175.0 * M_PI / 180.0;
        double theta_max_1 =  175.0 * M_PI / 180.0;

        vector<pair<double, double>> intervals_1;
        bool success_1 = solveMonotonicIntervals(
            A_n_1, B_n_1, C_n_1,
            A_d_1, B_d_1, C_d_1,
            a_t1, b_t1, c_t1,
            theta_min_1, theta_max_1,
            intervals_1
        );

        if (Delta1 > 0.0 && success_1) {
            // 循环型（图3a）
            printf("-------------theta1为循环型\n");

        } else if (Delta1 < 0.0) {
            // 单调型（图3b）
            printf("-------------theta1为单调型\n");
            for (auto &p : intervals_1) {
                printf("[%.4f, %.4f]\n", p.first, p.second);
            }
        } else {
            // 奇异型（图3c/d）
            printf("-------------theta1为奇异型\n");
        }

        //--------------求 theta2 可行范围----------
        double A_2 = -A_s(2,1);
        double B_2 = -B_s(2,1);
        double C_2 = -C_s(2,1);
        // 式37 38

        double sqrt_val_2 = sqrt(A_2 * A_2 + B_2 * B_2);

        // 注意：分母为 a = A_1，需判断是否为0
        double psi_minus_2, psi_max_2;

        psi_minus_2 = 2.0 * atan2(-B_2 - sqrt_val_2, A_2);   // 式(37)
        psi_max_2  = 2.0 * atan2(-B_2 + sqrt_val_2, A_2);   // 式(38)
        printf("-------------theta2为循环型\n");
        printf("theta2 对应比较最小值为：%f, 最大值为:%f \n",psi_minus_2,psi_max_2);


        //--------------求 theta3 可行范围----------

        double A_n_3 =  A_s(2,2);
        double B_n_3 =  B_s(2,2);
        double C_n_3 =  C_s(2,2);

        double A_d_3 = -A_s(2,0);
        double B_d_3 = -B_s(2,0);
        double C_d_3 = -C_s(2,0);

        double a_t3 = B_d_3 * C_n_3 - B_n_3 * C_d_3;
        double b_t3 = A_n_3 * C_d_3 - A_d_3 * C_n_3;
        double c_t3 = A_n_3 * B_d_3 - A_d_3 * B_n_3;

        double Delta3 = a_t3*a_t3 + b_t3*b_t3 - c_t3*c_t3;

        //确认theta3为循环（公式28 29） double check！！！！！！！！！！！！！
        double sqrtDelta3 = sqrt(Delta3);
        double psi_min3 = 2.0 * atan2(a_t3 - sqrtDelta3, b_t3 - c_t3);
        double psi_max3 = 2.0 * atan2(a_t3 + sqrtDelta3, b_t3 - c_t3);

        if (Delta3 > 0.0) {
            // 循环型（图3a）
            printf("-------------theta3为循环型\n");
            printf("theta3 对应比较最小值为：%f, 最大值为:%f \n",psi_min3,psi_max3);
        } else if (Delta3 < 0.0) {
            // 单调型（图3b）
            printf("-------------theta3为单调型\n");
        } else {
            // 奇异型（图3c/d）
            printf("-------------theta3为奇异型\n");
        }
    }



    // 4. 计算当前 arm_angle 下的 R_0_3
    Eigen::Matrix3d R_0_3 =
        A_s * std::sin(arm_angle)
        + B_s * std::cos(arm_angle)
        + C_s;

    // 5. 从 R_0_3 反解 theta_1, theta_2, theta_3
    //
    // 对于：
    // R_0_3 = Rz(theta_1) * Ry(-theta_2) * Rz(theta_3)
    //
    // 展开后有：
    // R(0,2) = -cos(theta_1) * sin(theta_2)
    // R(1,2) = -sin(theta_1) * sin(theta_2)
    // R(2,0) =  sin(theta_2) * cos(theta_3)
    // R(2,1) = -sin(theta_2) * sin(theta_3)
    // R(2,2) =  cos(theta_2)

    double sin_theta_2_abs = std::sqrt(
        R_0_3(2, 0) * R_0_3(2, 0)
        + R_0_3(2, 1) * R_0_3(2, 1)
    );

    // 主解：theta_2 in [0, pi]
    if (sin_theta_2_abs > eps) {
        theta_2 = std::atan2(sin_theta_2_abs, R_0_3(2, 2));

        theta_1 = std::atan2(
            -R_0_3(1, 2),
            -R_0_3(0, 2)
        );

        theta_3 = std::atan2(
            -R_0_3(2, 1),
            R_0_3(2, 0)
        );
    } else {
        // 奇异情况：sin(theta_2) 接近 0
        // 此时 theta_1 和 theta_3 耦合，无法唯一分开。
        // 这里保留 theta_3 = 0，把总的 z 方向旋转给 theta_1。
        theta_2 = std::atan2(0.0, R_0_3(2, 2));
        theta_3 = 0.0;

        if (R_0_3(2, 2) > 0.0) {
            // theta_2 ≈ 0，此时 R ≈ Rz(theta_1 + theta_3)
            theta_1 = std::atan2(R_0_3(1, 0), R_0_3(0, 0));
        } else {
            // theta_2 ≈ pi，此时也属于奇异，给一个可用分解
            theta_1 = std::atan2(-R_0_3(1, 0), -R_0_3(0, 0));
        }
    }

    theta_1 = normalizeAngle(theta_1);
    theta_2 = normalizeAngle(theta_2);
    theta_3 = normalizeAngle(theta_3);
    
    q_out[0] = (theta_1 + M_PI) * 180.0 / M_PI;// + M_PI 为offset
    q_out[1] = theta_2 * 180.0 / M_PI;
    q_out[2] = (theta_3 + M_PI) * 180.0 / M_PI;// + M_PI 为offset

    printf("theta_1 = %f deg\n", q_out[0] );
    printf("theta_2 = %f deg\n", q_out[1] );
    printf("theta_3 = %f deg\n", q_out[2] );
    printf("theta_4 = %f deg\n", q_out[3] );

    // 6. 验证一下分解误差
    Eigen::Matrix3d R_0_3_check =
        Rz(theta_1) * Ry_neg(theta_2) * Rz(theta_3);

    double R_0_3_err = (R_0_3_check - R_0_3).norm();


    // ================= 求 wrist joints: theta_5, theta_6, theta_7 =================
    //
    // 你的腕部轴定义：
    // R_4_7 = Rz(theta_5) * Ry(-theta_6) * Rz(theta_7)
    //
    // 其中第 6 关节是绕 y 负轴。

    double theta_5 = 0.0;
    double theta_6 = 0.0;
    double theta_7 = 0.0;

    {
        //-----------------------求解 Aw Bw Cw --------------------
        Eigen::Matrix3d A_w = R_3_4.transpose() * A_s.transpose() *R_0_desire;
        Eigen::Matrix3d B_w = R_3_4.transpose() * B_s.transpose() *R_0_desire;
        Eigen::Matrix3d C_w = R_3_4.transpose() * C_s.transpose() *R_0_desire;

        //---------------------求 theta5 的可行范围---------------------
        double A_n_5 = A_w(1,2);
        double B_n_5 = B_w(1,2);
        double C_n_5 = C_w(1,2);
        double A_d_5 = A_w(0,2);
        double B_d_5 = B_w(0,2);
        double C_d_5 = C_w(0,2);

        // 计算 a_t, b_t, c_t（公式27）
        double a_t5 = B_d_5 * C_n_5 - B_n_5 * C_d_5;
        double b_t5 = A_n_5 * C_d_5 - A_d_5 * C_n_5;
        double c_t5 = A_n_5 * B_d_5 - A_d_5 * B_n_5;

        // 判别式
        double Delta5 = a_t5*a_t5 + b_t5*b_t5 - c_t5*c_t5;

        //确认 theta5 为单调，按照关节范围对应最大最小
        double theta_min_5 = -175.0 * M_PI / 180.0;
        double theta_max_5 =  175.0 * M_PI / 180.0;

        vector<pair<double, double>> intervals_5;
        bool success_5 = solveMonotonicIntervals(
            A_n_5, B_n_5, C_n_5,
            A_d_5, B_d_5, C_d_5,
            a_t5,  b_t5,  c_t5,
            theta_min_5, theta_max_5,
            intervals_5
        );

        if (Delta5 > 0.0) {
            // 循环型（图3a）
            printf("-------------theta5为循环型\n");

        } else if (Delta5 < 0.0) {
            // 单调型（图3b）
            printf("-------------theta5为单调型\n");
            for (auto &p : intervals_5) {
                printf("[%.4f, %.4f]\n", p.first, p.second);
            }
        } else {
            // 奇异型（图3c/d）
            printf("-------------theta5为奇异型\n");
        }

        //-----------------------求 theta6 可行范围-------------------
            double A_6 = A_s(2,2);
            double B_6 = B_s(2,2);
            double C_6 = C_s(2,2);
            // 式37 38

            double sqrt_val_6 = sqrt(A_6 * A_6 + B_6 * B_6);

            // 注意：分母为 a = A_1，需判断是否为0
            double psi_minus_6, psi_max_6;

            psi_minus_6 = 2.0 * atan2(-B_6 - sqrt_val_6, A_6);   // 式(37)
            psi_max_6   = 2.0 * atan2(-B_6 + sqrt_val_6, A_6);   // 式(38)
            printf("-------------theta6为循环型\n");
            printf("theta6 对应比较最小值为：%f, 最大值为:%f \n",psi_minus_6,psi_max_6);


        //---------------------求 theta7 的可行范围---------------------
        double A_n_7 = A_w(2,1);
        double B_n_7 = B_w(2,1);
        double C_n_7 = C_w(2,1);
        double A_d_7 = -A_w(2,0);
        double B_d_7 = -B_w(2,0);
        double C_d_7 = -C_w(2,0);

        // 计算 a_t, b_t, c_t（公式27）
        double a_t7 = B_d_7 * C_n_7 - B_n_7 * C_d_7;
        double b_t7 = A_n_7 * C_d_7 - A_d_7 * C_n_7;
        double c_t7 = A_n_7 * B_d_7 - A_d_7 * B_n_7;

        // 判别式
        double Delta7 = a_t7*a_t7 + b_t7*b_t7 - c_t7*c_t7;

        //确认theta3为循环（公式28 29） double check！！！！！！！！！！！！！
        double sqrtDelta7 = sqrt(Delta7);
        double psi_min7 = 2.0 * atan2(a_t7 - sqrtDelta7, b_t7 - c_t7);
        double psi_max7 = 2.0 * atan2(a_t7 + sqrtDelta7, b_t7 - c_t7);

        if (Delta7 > 0.0) {
            // 循环型（图3a）
            printf("-------------theta7为循环型\n");
            printf("theta7 对应比较最小值为：%f, 最大值为:%f \n",psi_min7,psi_max7);
        } else if (Delta7 < 0.0) {
            // 单调型（图3b）
            printf("-------------theta7为单调型\n");
        } else {
            // 奇异型（图3c/d）
            printf("-------------theta7为奇异型\n");
        }
    }


    // R_0_desire 是目标末端旋转矩阵，也就是 ^0R_7^d
    // R_0_3 是你前面通过 arm_angle 算出来的 ^0R_3
    // R_3_4 是你由 theta_4 算出来的 ^3R_4
    Eigen::Matrix3d R_4_7 =
        R_3_4.transpose() * R_0_3.transpose() * R_0_desire;

    // 对于 R_4_7 = Rz(theta5) * Ry(-theta6) * Rz(theta7)
    //
    // R_4_7(0,2) = -cos(theta5) * sin(theta6)
    // R_4_7(1,2) = -sin(theta5) * sin(theta6)
    // R_4_7(2,0) =  sin(theta6) * cos(theta7)
    // R_4_7(2,1) = -sin(theta6) * sin(theta7)
    // R_4_7(2,2) =  cos(theta6)

    double sin_theta_6_abs = std::sqrt(
        R_4_7(2, 0) * R_4_7(2, 0)
        + R_4_7(2, 1) * R_4_7(2, 1)
    );

    if (sin_theta_6_abs > eps) {
        theta_6 = std::atan2(sin_theta_6_abs, R_4_7(2, 2));

        theta_5 = std::atan2(
            -R_4_7(1, 2),
            -R_4_7(0, 2)
        );

        theta_7 = std::atan2(
            -R_4_7(2, 1),
            R_4_7(2, 0)
        );
    } else {
        // 奇异情况：theta_6 接近 0 或 pi
        // 此时 theta_5 和 theta_7 耦合，不能唯一分开。
        // 常用处理：令 theta_7 = 0，把总旋转给 theta_5。
        theta_7 = 0.0;

        if (R_4_7(2, 2) > 0.0) {
            // theta_6 ≈ 0
            theta_6 = 0.0;

            // R ≈ Rz(theta_5 + theta_7)
            theta_5 = std::atan2(R_4_7(1, 0), R_4_7(0, 0));
        } else {
            // theta_6 ≈ pi
            theta_6 = M_PI;

            // 这里给一个可用分解
            theta_5 = std::atan2(-R_4_7(1, 0), -R_4_7(0, 0));
        }
    }

    theta_5 = normalizeAngle(theta_5);
    theta_6 = normalizeAngle(theta_6);
    theta_7 = normalizeAngle(theta_7);

    q_out[4] = (theta_5 + M_PI ) * 180.0 / M_PI;// + M_PI 为offset
    q_out[5] = theta_6 * 180.0 / M_PI;
    q_out[6] = (theta_7 + M_PI)* 180.0 / M_PI;// + M_PI 为offset

    printf("theta_5 = %f deg\n", q_out[4] );
    printf("theta_6 = %f deg\n", q_out[5] );
    printf("theta_7 = %f deg\n", q_out[6] );

    // 验证
    Eigen::Matrix3d R_4_7_check =
        Rz(theta_5) * Ry_neg(theta_6) * Rz(theta_7);

    double wrist_err = (R_4_7_check - R_4_7).norm();

    printf("R_0_3 decomposition error = %.12f\n", R_0_3_err);
    printf("R_4_7 decomposition error = %.12f\n", wrist_err);

}

/**
 * 将臂角限制到最近的可行区间
 * @param 当前臂角（rad）
 * @return best_psi 合理的臂角区间（rad）
 * 
 */
static inline double clampArmAngle(double psi)
{
    // 归一化到 [-π, π)
    while (psi < -M_PI) psi += 2.0 * M_PI;
    while (psi >=  M_PI) psi -= 2.0 * M_PI;

    // 已在区间1
    if (psi >= ARM_ANGLE_DISTRICT_1_1 && psi <= ARM_ANGLE_DISTRICT_1_2)
        return psi;

    // 已在区间2
    if (psi >= ARM_ANGLE_DISTRICT_2_1 && psi <= ARM_ANGLE_DISTRICT_2_2)
        return psi;

    // 不在任何区间 → 找最近的边界
    double dist_1_start = std::abs(psi - ARM_ANGLE_DISTRICT_1_1);
    double dist_1_end   = std::abs(psi - ARM_ANGLE_DISTRICT_1_2);
    double dist_2_start = std::abs(psi - ARM_ANGLE_DISTRICT_2_1);
    double dist_2_end   = std::abs(psi - ARM_ANGLE_DISTRICT_2_2);

    double min_dist = dist_1_start;
    double best_psi = ARM_ANGLE_DISTRICT_1_1;

    if (dist_1_end < min_dist) { min_dist = dist_1_end; best_psi = ARM_ANGLE_DISTRICT_1_2; }
    if (dist_2_start < min_dist) { min_dist = dist_2_start; best_psi = ARM_ANGLE_DISTRICT_2_1; }
    if (dist_2_end < min_dist) { min_dist = dist_2_end; best_psi = ARM_ANGLE_DISTRICT_2_2; }

    return best_psi;
}

// 只对关节2，4，6检查，其余关节权重为0
static double singularity_penalty(const double q[7]) {
    // 钟形惩罚：在奇异角度处峰值为1，sigma控制宽度
    // penalty = exp(-x^2 / (2*sigma^2))，sigma=0.15rad≈8.6度
    const double sigma2 = 0.15 * 0.15;
    const int sing_idx[] = {1, 3, 5};
    const double sing_angles[] = {0.0, M_PI, -M_PI};  // 每个关节的奇异角

    double penalty = 0.0;
    for (int k = 0; k < 3; k++) {
        int j = sing_idx[k];
        for (int a = 0; a < 3; a++) {
            double diff = q[j] - sing_angles[a];
            penalty += exp(-diff * diff / (2.0 * sigma2));
        }
    }
    return penalty;
}




/**
 * @param q 计算后的角度
 * @param q_init 当前初始角度
 * @return score 打分系统
 */
static double score_solution(const double q[7],const double q_init[7],
                              const double q_prev[7]) {
    const double w_dq   = 1.0;
    const double w_vel  = 0.5;
    const double w_sing = 5.0;

    double score = 0.0;

    for (int i = 0; i < 7; i++) {
        // 超出限位直接淘汰
        if (q[i] < Q_MIN[i]*180/M_PI || q[i] > Q_MAX[i]*180/M_PI) return INFINITY;

        // 1. 关节角度变化
        double dq = q[i] - q_init[i];
        score += w_dq * dq * dq;

        // 2. 速度变化（加速度代理）：(v_new - v_old) 其中 v = dq/dt
        double v_new = (q[i]    - q_init[i]) / DT;  // 当前步速度
        double v_old = (q_init[i] - q_prev[i]) / DT; // 上一步速度
        double dv = v_new - v_old;
        score += w_vel * dv * dv;}

    // 3. 奇异性（钟形，只在奇异角度附近才大）
    score += w_sing * singularity_penalty(q);

    return score;
}

// 主函数：给定目标位姿，计算臂角可行区间
// 返回可行区间列表（已是各关节约束的交集）
vector<pair<double,double>> compute_arm_angle_intervals(const double T_target[4][4])
{
    Eigen::Matrix3d R_0_d;
    Eigen::Vector3d P_0_d;
    getPoseFromArray(T_target, R_0_d, P_0_d);

    Eigen::Vector3d l_bs(0, 0, dh[0][2]);
    Eigen::Vector3d l_wt(0, 0, sqrt(dh[6][2]*dh[6][2] + dh[5][2]*dh[5][2]));
    Eigen::Vector3d x_sw = P_0_d - l_bs - R_0_d * l_wt;

    double norm_sw = x_sw.norm();
    Eigen::Vector3d u_sw = Eigen::Vector3d::Zero();
    if (norm_sw > 1e-12) u_sw = x_sw / norm_sw;

    // theta4（与psi无关）
    double cos4 = (norm_sw*norm_sw - dh[2][2]*dh[2][2] - dh[4][2]*dh[4][2])
                  / (2.0 * dh[2][2] * dh[4][2]);
    cos4 = std::clamp(cos4, -1.0, 1.0);
    double theta_4 = -acos(cos4); // 默认取负分支，与主IK一致
    printf("theta_4 = %.4f deg\n", theta_4*180/M_PI);

    // R_3_4
    double c4 = cos(theta_4), s4 = sin(theta_4);
    Eigen::Matrix3d R_3_4;
    R_3_4 << c4, 0, -s4,  0, 1, 0,  s4, 0, c4;

    // 求 theta_1_ref, theta_2_ref（psi=0 参考位形）
    Eigen::Vector3d l_se(0, 0, dh[2][2]);
    Eigen::Vector3d l_ew(0, 0, dh[4][2]);
    Eigen::Vector3d p_ref = l_se + R_3_4 * l_ew; // R_2_3=I

    double px = p_ref(0), py = p_ref(1), pz = p_ref(2);
    double x = x_sw(0), y = x_sw(1), z = x_sw(2);
    double r = sqrt(px*px + pz*pz);

    double theta_1_ref = 0, theta_2_ref = 0;
    if (r > 1e-12) {
        double alpha = atan2(px, pz);
        double beta  = acos(std::clamp(z/r, -1.0, 1.0));

        auto try_cand = [&](double t2) -> pair<double,double> {
            double c2 = cos(t2), s2 = sin(t2);
            double a = c2*px - s2*pz;
            double t1 = atan2(y, x) - atan2(py, a);
            while (t1 > M_PI) t1 -= 2*M_PI;
            while (t1 <-M_PI) t1 += 2*M_PI;
            return {t1, t2};
        };

        auto [t1a, t2a] = try_cand(alpha + beta);
        auto [t1b, t2b] = try_cand(alpha - beta);

        // 用FK残差选
        auto Rz = [](double th) -> Eigen::Matrix3d {
            Eigen::Matrix3d R; double c=cos(th),s=sin(th);
            R << c,-s,0, s,c,0, 0,0,1; return R;
        };
        auto Ry_neg = [](double th) -> Eigen::Matrix3d {
            Eigen::Matrix3d R; double c=cos(th),s=sin(th);
            R << c,0,-s, 0,1,0, s,0,c; return R;
        };
        double ea = (Rz(t1a)*Ry_neg(t2a)*p_ref - x_sw).norm();
        double eb = (Rz(t1b)*Ry_neg(t2b)*p_ref - x_sw).norm();
        if (ea <= eb) { theta_1_ref = t1a; theta_2_ref = t2a; }
        else          { theta_1_ref = t1b; theta_2_ref = t2b; }
    }

    // 构造 A_s, B_s, C_s（肩部）
    Eigen::Matrix3d u_cross;
    u_cross << 0, -u_sw(2), u_sw(1),
               u_sw(2), 0, -u_sw(0),
              -u_sw(1), u_sw(0), 0;
    auto Rz_l = [](double th) -> Eigen::Matrix3d {
        Eigen::Matrix3d R; double c=cos(th),s=sin(th);
        R << c,-s,0, s,c,0, 0,0,1; return R;
    };
    auto Ry_neg_l = [](double th) -> Eigen::Matrix3d {
        Eigen::Matrix3d R; double c=cos(th),s=sin(th);
        R << c,0,-s, 0,1,0, s,0,c; return R;
    };
    Eigen::Matrix3d R_0_3_ref = Rz_l(theta_1_ref) * Ry_neg_l(theta_2_ref);
    Eigen::Matrix3d A_s = u_cross * R_0_3_ref;
    Eigen::Matrix3d B_s = -u_cross * u_cross * R_0_3_ref;
    Eigen::Matrix3d C_s = (u_sw * u_sw.transpose()) * R_0_3_ref;

    // 腕部矩阵
    Eigen::Matrix3d A_w = R_3_4.transpose() * A_s.transpose() * R_0_d;
    Eigen::Matrix3d B_w = R_3_4.transpose() * B_s.transpose() * R_0_d;
    Eigen::Matrix3d C_w = R_3_4.transpose() * C_s.transpose() * R_0_d;

    const double th_min = Q_MIN[0]; // ±175° in rad，所有关节相同
    const double th_max = Q_MAX[0];

    // ---- theta1: atan2(-R_0_3(1,2), -R_0_3(0,2)) → atan2型 ----
    auto iv1 = atan2_joint_interval(
        -A_s(1,2), -B_s(1,2), -C_s(1,2),
        -A_s(0,2), -B_s(0,2), -C_s(0,2),
        th_min, th_max);

    // ---- theta2: acos(R_0_3(2,2)) = acos(A_s(2,2)*sin+B_s(2,2)*cos+C_s(2,2)) ----
    // R_0_3(2,2) = cos(theta2)，theta2 in [0,pi] → theta_min=0, theta_max=pi 与关节限位取min
    auto iv2 = acos_joint_interval(
        B_s(2,2), A_s(2,2), C_s(2,2), // cos型：A*cos+B*sin+C → swap A,B for acos helper
        std::max(th_min, 0.0), std::min(th_max, M_PI));

    // ---- theta3: atan2(-R(2,1)/s2, R(2,0)/s2) → atan2型 ----
    auto iv3 = atan2_joint_interval(
        -A_s(2,1), -B_s(2,1), -C_s(2,1),
         A_s(2,0),  B_s(2,0),  C_s(2,0),
        th_min, th_max);

    // ---- theta5: atan2(-R_4_7(1,2)/-R_4_7(0,2)) → atan2型 ----
    auto iv5 = atan2_joint_interval(
        -A_w(1,2), -B_w(1,2), -C_w(1,2),
        -A_w(0,2), -B_w(0,2), -C_w(0,2),
        th_min, th_max);

    // ---- theta6: acos(R_4_7(2,2)) ----
    auto iv6 = acos_joint_interval(
        B_w(2,2), A_w(2,2), C_w(2,2),
        std::max(th_min, 0.0), std::min(th_max, M_PI));

    // ---- theta7: atan2(-R_4_7(2,1), R_4_7(2,0)) → atan2型 ----
    auto iv7 = atan2_joint_interval(
        -A_w(2,1), -B_w(2,1), -C_w(2,1),
         A_w(2,0),  B_w(2,0),  C_w(2,0),
        th_min, th_max);

    auto print_iv = [](const char* name, const char* type,
                       const vector<pair<double,double>>& iv) {
        printf("  %-8s [%s] %zu segment(s):", name, type, iv.size());
        if (iv.empty()) printf("  EMPTY (joint unreachable)");
        for (auto& s : iv)
            printf("  [%7.2f, %7.2f] deg", s.first*180/M_PI, s.second*180/M_PI);
        printf("\n");
    };
    printf("\n===== 各关节ψ可行区间 =====\n");
    print_iv("theta1", "atan2", iv1);
    print_iv("theta2", "acos ", iv2);
    print_iv("theta3", "atan2", iv3);
    print_iv("theta4", "const", {});  // 常数，不参与约束
    print_iv("theta5", "atan2", iv5);
    print_iv("theta6", "acos ", iv6);
    print_iv("theta7", "atan2", iv7);

    // 求所有关节约束的交集
    auto res = full_range();
    for (auto* iv : {&iv1, &iv2, &iv3, &iv5, &iv6, &iv7})
        res = intersect(res, *iv);

    return res;
}


/**
 * @param T_target 目标位姿
 * @param q_init 当前角度（degree）
 * @param psi 当前臂角（rad）
 * @param q_out 角度（degree）
 * 
 */
void analytical_ik_paper(const double T_target[4][4], const double q_init[7], double psi, double q_out[7]){
    double arm_angle = psi;//37.904784*M_PI/180

    Eigen::Matrix3d R_0_desire;
    Eigen::Vector3d P_0_desire;
    //获取base下的SW坐标,及其单位向量
    Eigen::Vector3d x_0_sw; // 在 base frame中 获取sw向量
    Eigen::Vector3d u_0_sw;//x_sw_0 的单位向量
    double q_input[7];
    memcpy(q_input,q_init,7*sizeof(double));

    getPoseFromArray(T_target,R_0_desire,P_0_desire);//获取目标位姿的posi rot
    Eigen::Vector3d l_0_bs(0,0,dh[0][2]);
    Eigen::Vector3d l_7_wt(0,0,sqrt(pow(dh[6][2],2) + pow(dh[5][2],2)));
    x_0_sw = P_0_desire - l_0_bs - R_0_desire *l_7_wt;
    printf("x_0_sw的位置为: %f,%f,%f\n",x_0_sw(0),x_0_sw(1),x_0_sw(2));

    //在此先检验一遍 奇异
    sigualrity_type is_singular = check_near_singularity(q_input,x_0_sw);
    if(is_singular<0){
        printf("------------------------在singular附近------------------------- \n");
    }
    
    double norm_x_0_sw = x_0_sw.norm();
    if (norm_x_0_sw > 1e-12) {
        u_0_sw = x_0_sw / norm_x_0_sw;
    } else {
        u_0_sw.setZero();  
    }

    //==============解 theta 4=========== checked
    double cos_theta4 = (pow(norm_x_0_sw,2) - pow(dh[2][2],2) - pow(dh[4][2],2)) / ( 2* dh[2][2] * dh[4][2]);
    
    if (cos_theta4 < -1.0) cos_theta4 = -1.0;
    if (cos_theta4 > 1.0) cos_theta4 = 1.0;
    double acos_val = std::acos(cos_theta4);
    // q_out[3] = -theta_4 * 180/pi，所以 theta_4=-acos → q_out[3]=+|v|，theta_4=+acos → q_out[3]=-|v|
    double q4_cand_pos =  acos_val * 180.0 / M_PI;  // theta_4 = -acos → q_out[3] 为正
    double q4_cand_neg = -acos_val * 180.0 / M_PI;  // theta_4 = +acos → q_out[3] 为负
    double theta_4 = (std::abs(unwrap_delta(q4_cand_pos, q_input[3])) <= std::abs(unwrap_delta(q4_cand_neg, q_input[3])))
                     ? -acos_val : acos_val;
    printf("------------第四个关节角度为：%f \n",theta_4*180/M_PI);
    q_out[3] = - theta_4 * 180.0 / M_PI;//负号为offset

    // ===== 构造 R_3_4 (theta_4 已求出，提到前面供后面复用) =====
    double c4 = std::cos(theta_4);
    double s4 = std::sin(theta_4);

    //暂时不使用fk计算（坐标系不一样），由于 j4是绕着 base frame 的 -y轴旋转的，按照旋转为：
    Eigen::Matrix3d R_3_4;
    R_3_4 <<  c4, 0.0, -s4,
            0.0, 1.0, 0.0,
            s4, 0.0,  c4;


    auto normalizeAngle = [](double angle) {
        while (angle > M_PI) angle -= 2.0 * M_PI;
        while (angle < -M_PI) angle += 2.0 * M_PI;
        return angle;
    };
    //第1个轴是绕base frame 的 z
    auto Rz = [](double theta) {
        double c = std::cos(theta);
        double s = std::sin(theta);
        Eigen::Matrix3d R;
        R << c, -s, 0.0,
             s,  c, 0.0,
             0.0, 0.0, 1.0;
        return R;
    };
    //第二个轴是绕base frame 的 -y
    auto Ry_neg = [](double theta) {
        double c = std::cos(theta);
        double s = std::sin(theta);
        Eigen::Matrix3d R;
        R <<  c, 0.0, -s,
             0.0, 1.0, 0.0,
             s, 0.0,  c;
        return R;
    };
    const double eps = 1e-9;

    //求当臂角为0时候的theta_1_ref plane 和theta_1_ref plane 对应公式 14
    double theta_1_ref = 0.0;
    double theta_2_ref = 0.0;

    // ================= 求 theta_1_ref, theta_2_ref =================
    //
    // 目标形式：
    // x_0_sw = Rz(theta_1_ref) * Ry(-theta_2_ref) * p
    //
    // p = R_2_3(theta3 = 0) * (l_3_se + R_3_4(theta4) * l_4_ew)
    Eigen::Matrix3d R_2_3_theta3_zero = Eigen::Matrix3d::Identity();
    Eigen::Vector3d l_3_se(0.0, 0.0, dh[2][2]);
    Eigen::Vector3d l_4_ew(0.0, 0.0, dh[4][2]);

    // 已知向量 p
    Eigen::Vector3d p = R_2_3_theta3_zero * (l_3_se + R_3_4 * l_4_ew);

    // 现在解：x_0_sw = Rz(theta_1_ref) * Ry(-theta_2_ref) * p
    // [ c1·c2,  -s1,  -c1·s2 ]                                   
    // [ s1·c2,   c1,  -s1·s2 ]                                   
    // [ s2,       0,     c2   ]   
    //有两个未知数，3组方程。直接解耦求解

    const double eps_1 = 1e-12;

    const double x = x_0_sw(0);
    const double y = x_0_sw(1);
    const double z = x_0_sw(2);

    const double px = p(0);
    const double py = p(1);
    const double pz = p(2);

    const double r = std::sqrt(px * px + pz * pz);

    if (r < eps_1) {
        printf("求 theta_2_ref 失败：p 在 y 轴附近，出现奇异。\n");
    } else {
        double cos_arg = z / r;

        if (cos_arg > 1.0) cos_arg = 1.0;
        if (cos_arg < -1.0) cos_arg = -1.0;

        const double alpha = std::atan2(px, pz);
        const double beta = std::acos(cos_arg);

        double theta_2_candidate_1 = alpha + beta;
        double theta_2_candidate_2 = alpha - beta;

        double theta_1_candidate_1 = 0.0;
        double theta_1_candidate_2 = 0.0;

        {
            double c2 = std::cos(theta_2_candidate_1);
            double s2 = std::sin(theta_2_candidate_1);

            double a = c2 * px - s2 * pz;

            theta_1_candidate_1 =
                std::atan2(y, x) - std::atan2(py, a);
        }

        {
            double c2 = std::cos(theta_2_candidate_2);
            double s2 = std::sin(theta_2_candidate_2);

            double a = c2 * px - s2 * pz;

            theta_1_candidate_2 =
                std::atan2(y, x) - std::atan2(py, a);
        }

        theta_1_candidate_1 = normalizeAngle(theta_1_candidate_1);
        theta_2_candidate_1 = normalizeAngle(theta_2_candidate_1);

        theta_1_candidate_2 = normalizeAngle(theta_1_candidate_2);
        theta_2_candidate_2 = normalizeAngle(theta_2_candidate_2);

        // 用正运动学残差选一个更匹配的解
        Eigen::Vector3d x_check_1 =
            Rz(theta_1_candidate_1) * Ry_neg(theta_2_candidate_1) * p;

        Eigen::Vector3d x_check_2 =
            Rz(theta_1_candidate_2) * Ry_neg(theta_2_candidate_2) * p;

        double err_1 = (x_check_1 - x_0_sw).norm();
        double err_2 = (x_check_2 - x_0_sw).norm();

        if (err_1 <= err_2) {
            theta_1_ref = theta_1_candidate_1;
            theta_2_ref = theta_2_candidate_1;
        } else {
            theta_1_ref = theta_1_candidate_2;
            theta_2_ref = theta_2_candidate_2;
        }

        printf("theta_1_ref candidate 1: %f deg, theta_2_ref candidate 1: %f deg, err: %.12f\n",
            theta_1_candidate_1 * 180.0 / M_PI,
            theta_2_candidate_1 * 180.0 / M_PI,
            err_1);

        printf("theta_1_ref candidate 2: %f deg, theta_2_ref candidate 2: %f deg, err: %.12f\n",
            theta_1_candidate_2 * 180.0 / M_PI,
            theta_2_candidate_2 * 180.0 / M_PI,
            err_2);

        printf("selected theta_1_ref: %f deg\n", theta_1_ref * 180.0 / M_PI);
        printf("selected theta_2_ref: %f deg\n", theta_2_ref * 180.0 / M_PI);
    }
    // ================= 根据 arm_angle 求 theta_1, theta_2, theta_3 =================
    //
    // 这里假设肩部前三个旋转为：
    // R_0_1 = Rz(theta_1)
    // R_1_2 = Ry(-theta_2)
    // R_2_3 = Rz(theta_3)
    //
    // 所以：
    // R_0_3 = Rz(theta_1) * Ry(-theta_2) * Rz(theta_3)

    double theta_1 = 0.0;
    double theta_2 = 0.0;
    double theta_3 = 0.0;

    // 1. 先构造 arm_angle = 0 时的参考 R_0_3_ref
    // theta_3_ref = 0
    Eigen::Matrix3d R_0_3_ref =
        Rz(theta_1_ref) * Ry_neg(theta_2_ref) * Rz(0.0);

    // 2. 构造 [u_0_sw x] 反对称矩阵
    Eigen::Matrix3d u_cross;
    u_cross << 0.0,        -u_0_sw(2),  u_0_sw(1),
            u_0_sw(2),   0.0,       -u_0_sw(0),
            -u_0_sw(1),   u_0_sw(0),  0.0;

    // 3. 根据论文公式构造 A_s, B_s, C_s
    Eigen::Matrix3d A_s = u_cross * R_0_3_ref;
    Eigen::Matrix3d B_s = -u_cross * u_cross * R_0_3_ref;
    Eigen::Matrix3d C_s = (u_0_sw * u_0_sw.transpose()) * R_0_3_ref;

    // 腕部矩阵
    Eigen::Matrix3d A_w = R_3_4.transpose() * A_s.transpose() * R_0_desire;
    Eigen::Matrix3d B_w = R_3_4.transpose() * B_s.transpose() * R_0_desire;
    Eigen::Matrix3d C_w = R_3_4.transpose() * C_s.transpose() * R_0_desire;

    // 4. 计算当前 arm_angle 下的 R_0_3
    Eigen::Matrix3d R_0_3 =
        A_s * std::sin(arm_angle)
        + B_s * std::cos(arm_angle)
        + C_s;
    

    // R_0_desire 是目标末端旋转矩阵，也就是 ^0R_7^d
    // R_0_3 是你前面通过 arm_angle 算出来的 ^0R_3
    // R_3_4 是你由 theta_4 算出来的 ^3R_4
    Eigen::Matrix3d R_4_7 =
        R_3_4.transpose() * R_0_3.transpose() * R_0_desire;

    // 5. 从 R_0_3 反解 theta_1, theta_2, theta_3
    //
    // 对于：
    // R_0_3 = Rz(theta_1) * Ry(-theta_2) * Rz(theta_3)
    //
    // 展开后有：
    // R(0,2) = -cos(theta_1) * sin(theta_2)
    // R(1,2) = -sin(theta_1) * sin(theta_2)
    // R(2,0) =  sin(theta_2) * cos(theta_3)
    // R(2,1) = -sin(theta_2) * sin(theta_3)
    // R(2,2) =  cos(theta_2)

    double sin_theta_2_abs = std::sqrt(
        R_0_3(2, 0) * R_0_3(2, 0)
        + R_0_3(2, 1) * R_0_3(2, 1)
    );

    // 主解：theta_2 in [0, pi]
    if (sin_theta_2_abs > eps) {
        theta_2 = std::atan2(sin_theta_2_abs, R_0_3(2, 2));

        theta_1 = std::atan2(
            -R_0_3(1, 2),
            -R_0_3(0, 2)
        );

        theta_3 = std::atan2(
            -R_0_3(2, 1),
            R_0_3(2, 0)
        );
    } else {
        // 奇异情况：sin(theta_2) 接近 0
        // 此时 theta_1 和 theta_3 耦合，无法唯一分开。
        // 这里保留 theta_3 = 0，把总的 z 方向旋转给 theta_1。
        theta_2 = std::atan2(0.0, R_0_3(2, 2));
        theta_3 = 0.0;

        if (R_0_3(2, 2) > 0.0) {
            // theta_2 ≈ 0，此时 R ≈ Rz(theta_1 + theta_3)
            theta_1 = std::atan2(R_0_3(1, 0), R_0_3(0, 0));
        } else {
            // theta_2 ≈ pi，此时也属于奇异，给一个可用分解
            theta_1 = std::atan2(-R_0_3(1, 0), -R_0_3(0, 0));
        }
    }

    theta_1 = normalizeAngle(theta_1);
    theta_2 = normalizeAngle(theta_2);
    theta_3 = normalizeAngle(theta_3);
    
    q_out[0] = (theta_1 + M_PI) * 180.0 / M_PI;// + M_PI 为offset
    q_out[1] = theta_2 * 180.0 / M_PI;
    q_out[2] = (theta_3 + M_PI) * 180.0 / M_PI;// + M_PI 为offset

    printf("theta_1 = %f deg\n", q_out[0] );
    printf("theta_2 = %f deg\n", q_out[1] );
    printf("theta_3 = %f deg\n", q_out[2] );
    printf("theta_4 = %f deg\n", q_out[3] );

    // 6. 验证一下分解误差
    Eigen::Matrix3d R_0_3_check =
        Rz(theta_1) * Ry_neg(theta_2) * Rz(theta_3);

    double R_0_3_err = (R_0_3_check - R_0_3).norm();


    // ================= 求 wrist joints: theta_5, theta_6, theta_7 =================
    //
    // 你的腕部轴定义：
    // R_4_7 = Rz(theta_5) * Ry(-theta_6) * Rz(theta_7)
    //
    // 其中第 6 关节是绕 y 负轴。

    double theta_5 = 0.0;
    double theta_6 = 0.0;
    double theta_7 = 0.0;

    
    // 对于 R_4_7 = Rz(theta5) * Ry(-theta6) * Rz(theta7)
    //
    // R_4_7(0,2) = -cos(theta5) * sin(theta6)
    // R_4_7(1,2) = -sin(theta5) * sin(theta6)
    // R_4_7(2,0) =  sin(theta6) * cos(theta7)
    // R_4_7(2,1) = -sin(theta6) * sin(theta7)
    // R_4_7(2,2) =  cos(theta6)

    double sin_theta_6_abs = std::sqrt(
        R_4_7(2, 0) * R_4_7(2, 0)
        + R_4_7(2, 1) * R_4_7(2, 1)
    );

    if (sin_theta_6_abs > eps) {
        theta_6 = std::atan2(sin_theta_6_abs, R_4_7(2, 2));

        theta_5 = std::atan2(
            -R_4_7(1, 2),
            -R_4_7(0, 2)
        );

        theta_7 = std::atan2(
            -R_4_7(2, 1),
            R_4_7(2, 0)
        );
    } else {
        // 奇异情况：theta_6 接近 0 或 pi
        // 此时 theta_5 和 theta_7 耦合，不能唯一分开。
        // 常用处理：令 theta_7 = 0，把总旋转给 theta_5。
        theta_7 = 0.0;

        if (R_4_7(2, 2) > 0.0) {
            // theta_6 ≈ 0
            theta_6 = 0.0;

            // R ≈ Rz(theta_5 + theta_7)
            theta_5 = std::atan2(R_4_7(1, 0), R_4_7(0, 0));
        } else {
            // theta_6 ≈ pi
            theta_6 = M_PI;

            // 这里给一个可用分解
            theta_5 = std::atan2(-R_4_7(1, 0), -R_4_7(0, 0));
        }
    }

    theta_5 = normalizeAngle(theta_5);
    theta_6 = normalizeAngle(theta_6);
    theta_7 = normalizeAngle(theta_7);

    q_out[4] = (theta_5 + M_PI ) * 180.0 / M_PI;// + M_PI 为offset
    q_out[5] = theta_6 * 180.0 / M_PI;
    q_out[6] = (theta_7 + M_PI)* 180.0 / M_PI;// + M_PI 为offset

    printf("theta_5 = %f deg\n", q_out[4] );
    printf("theta_6 = %f deg\n", q_out[5] );
    printf("theta_7 = %f deg\n", q_out[6] );

    // 验证
    Eigen::Matrix3d R_4_7_check =
        Rz(theta_5) * Ry_neg(theta_6) * Rz(theta_7);

    double wrist_err = (R_4_7_check - R_4_7).norm();

    printf("R_0_3 decomposition error = %.12f\n", R_0_3_err);
    printf("R_4_7 decomposition error = %.12f\n", wrist_err);

}


static double eval_psi(double psi,
                       const double T_target[4][4],
                       const double q_init[N_JOINTS],
                       const double q_prev[N_JOINTS],
                       double q_out[N_JOINTS]) {
    analytical_ik_paper(T_target, q_init, psi, q_out);
    for (int j = 0; j < N_JOINTS; j++)
        if (!std::isfinite(q_out[j])) return INFINITY;
    return score_solution(q_out, q_init, q_prev);
}

// 黄金分割 1D 搜索，在 [lo, hi] 内找最小值
// 返回最低得分，psi_best 为最优臂角，q_best 为对应关节角
static double golden_section_1d(double lo, double hi, double tol,
                                const double T_target[4][4],
                                const double q_init[N_JOINTS],
                                const double q_prev[N_JOINTS],
                                double q_best[N_JOINTS],
                                double& psi_best) {
    const double phi   = 0.6180339887498949;   // 1/φ
    const double phi_c = 1.0 - phi;             // ~0.382

    double a = lo, b = hi;
    double x1 = a + phi_c * (b - a);   // 左内点
    double x2 = a + phi   * (b - a);   // 右内点

    double q1[N_JOINTS], q2[N_JOINTS];
    double f1 = eval_psi(x1, T_target, q_init, q_prev, q1);
    double f2 = eval_psi(x2, T_target, q_init, q_prev, q2);
    int evals = 2;

    while ((b - a) > tol && evals < 20) {
        if (f1 < f2) {
            // 抛弃右段 [x2, b]，旧 x1 变成新 x2（复用）
            b  = x2;
            x2 = x1;   f2 = f1;
            x1 = a + phi_c * (b - a);
            f1 = eval_psi(x1, T_target, q_init, q_prev, q1);
        } else {
            // 抛弃左段 [a, x1]，旧 x2 变成新 x1（复用）
            a  = x1;
            x1 = x2;   f1 = f2;
            x2 = a + phi * (b - a);
            f2 = eval_psi(x2, T_target, q_init, q_prev, q2);
        }
        evals++;
    }

    psi_best = (a + b) / 2.0;
    return eval_psi(psi_best, T_target, q_init, q_prev, q_best);
}

// 对外接口：与 select_optimal_ik 相同签名
int select_optimal_ik_golden(const double T_target[4][4],
                             const double q_init[N_JOINTS],
                             const double q_prev[N_JOINTS],
                             double q_best[N_JOINTS]) {
    const double search_radius = 5.0 * (M_PI / 180.0);   // ±5°
    const double tol = 0.005 * (M_PI / 180.0);            // ~0.005°

    double psi_center = arm_plane_angle(q_init);
    printf("[黄金分割] 当前臂角: %.3f°\n", psi_center * 180.0 / M_PI);

    double best_score = INFINITY;
    int found = 0;

    const double iv[2][2] = {
        {ARM_ANGLE_DISTRICT_1_1, ARM_ANGLE_DISTRICT_1_2},
        {ARM_ANGLE_DISTRICT_2_1, ARM_ANGLE_DISTRICT_2_2}
    };

    for (int k = 0; k < 2; k++) {
        double lo = std::max(iv[k][0], psi_center - search_radius);
        double hi = std::min(iv[k][1], psi_center + search_radius);
        if (hi - lo < tol) continue;

        printf("[黄金分割] 区间%d [%.2f°, %.2f°] 搜索中...\n",
               k, lo * 180.0 / M_PI, hi * 180.0 / M_PI);

        double psi_best_k;
        double q_k[N_JOINTS];
        double score_k = golden_section_1d(lo, hi, tol,
                                           T_target, q_init, q_prev,
                                           q_k, psi_best_k);

        printf("[黄金分割] 区间%d 最优臂角: %.4f°, 得分: %.3f\n",
               k, psi_best_k * 180.0 / M_PI, score_k);

        if (score_k < best_score) {
            best_score = score_k;
            memcpy(q_best, q_k, N_JOINTS * sizeof(double));
            found = 1;
        }
    }

    return found;
}




// q_prev: 上上时刻关节角（用于估计加速度）
/**
 * @param q_init 上一时刻（t-1）的关节角
 * @param q_prev 上上时刻（t-2）的关节角
 * @param psi_center 上一时刻（t-1）的臂角
 * @param q_best 当前时刻（t）的候选解
 */
int select_optimal_ik(const double T_target[4][4],
                      const double q_init[N_JOINTS],
                      const double q_prev[N_JOINTS],
                      double q_best[N_JOINTS]) {
    const double step = 0.6 * (M_PI / 180.0);

    double best_score = INFINITY;
    int found = 0;

    double psi_center = arm_plane_angle(q_init);
    printf("*******************理想采样臂角为：%f \n",psi_center*180/M_PI);

    for (int i = 0; i < 15; i++) {
        double psi_sample = psi_center + (-7 + i) * step;  // 7点对称采样
        double psi = clampArmAngle(psi_sample);//限制臂角
        printf("*******************采样臂角为：%f \n",psi*180/M_PI);
        double q_cand[N_JOINTS];
        analytical_ik_paper(T_target, q_init, psi, q_cand);

        // NaN/Inf 检查
        int valid = 1;
        for (int j = 0; j < N_JOINTS; j++)
            if (!isfinite(q_cand[j])) { valid = 0; break; }
        if (!valid) continue;

        double s = score_solution(q_cand, q_init, q_prev);
        if (s < best_score) {
            best_score = s;
            memcpy(q_best, q_cand, N_JOINTS * sizeof(double));
            found = 1;
        }
    }
    

    return found;
}






int main(){
    // analytical_fk();/
    double T_target[4][4] = 
        {{-0.02385176519, 0.7667205121, -0.6415378006, 251.5004571},
        {0.05935864847, 0.6416743366, 0.7646767922, 885.2382746},
        {0.9979517244, -0.01984192549, -0.06081655896, 316.4981235},
        {0, 0, 0, 1}};
    auto start = std::chrono::high_resolution_clock::now();
    double joint_angle[7] = {};
    double joint_init[7] = {52, 65, 42, 51, 62, 50.5, 85};
    double joint_outcome[7] = {50, 60, 40, 50, 60, 50, 80};

    

    // select_optimal_ik_golden(T_target,joint_init,joint_init,joint_angle);

    double cur_arm_angel = arm_plane_angle(joint_outcome);
    printf("************当前臂角为： %f\n",cur_arm_angel);
    analytical_ik_paper_with_arm_angle_cal(T_target,joint_angle);

    for(int i = 0;i<7;i++){
        printf("joint_angle%d的角度为：%f\n",i,joint_angle[i]);
    }
    double T_compare[4][4] = {};
    forward_kinematics(joint_angle,T_compare,7);
    compareTransforms(T_target,T_compare);
    for (int i = 0; i < 4; ++i) {
        std::cout << "  ";
        for (int j = 0; j < 4; ++j) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(4) << T_compare[i][j];
        }
        std::cout << "\n";
    }
    std::cout << std::endl;


    printf("综上所述，存在的臂角区间为：[ %f, %f ]U[ %f, %f ]°\n ",
        ARM_ANGLE_DISTRICT_1_1*180/M_PI,ARM_ANGLE_DISTRICT_1_2*180/M_PI, 
        ARM_ANGLE_DISTRICT_2_1*180/M_PI,ARM_ANGLE_DISTRICT_2_2*180/M_PI);


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> elapsed = end - start;
    printf("\n===== 用时: %.3f us (%.3f ms) =====\n", elapsed.count(), elapsed.count() / 1000.0);

    return 0;
}
//编译： g++ analytical_ik.cpp r7_fk_lib.cpp -o analytical_ik

//综上所述，存在的臂角区间为：【0.000000  0.9021】U[1.1019   2.646297] 
