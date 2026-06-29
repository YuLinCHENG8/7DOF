#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <Eigen/Geometry>
#include <iomanip>

using std::vector;
using std::pair;
using namespace std;

// 关节限位（弧度）
static const double Q_MIN[7] = {-3.054326,-3.054326,-3.054326,-3.054326,-3.054326,-3.054326,-3.054326};
static const double Q_MAX[7] = { 3.054326, 3.054326, 3.054326, 3.054326, 3.054326, 3.054326, 3.054326};

const double dh[7][4] = {
    {  0, 0, 84+95,   0},
    {-90, 0,     0,   0},
    { 90, 0, 215+260, 0},
    {-90, 0,     0,   0},
    { 90, 0, 415+60,  0},
    {-90, 0,  -2.5,   0},
    { 90, 0,   145,   0},
};

void forward_kinematics(const double joint_angles[7], double T_out[4][4], int joint = 7);

void getPoseFromArray(const double T[4][4], Eigen::Matrix3d& R, Eigen::Vector3d& p) {
    R << T[0][0], T[0][1], T[0][2],
         T[1][0], T[1][1], T[1][2],
         T[2][0], T[2][1], T[2][2];
    p << T[0][3], T[1][3], T[2][3];
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

int main() {
    double T_target[4][4] = {
        {-0.02385176519, 0.7667205121, -0.6415378006, 251.5004571},
        { 0.05935864847, 0.6416743366,  0.7646767922, 885.2382746},
        { 0.9979517244, -0.01984192549,-0.06081655896, 316.4981235},
        {0, 0, 0, 1}
    };

    auto intervals = compute_arm_angle_intervals(T_target);

    printf("\n===== 臂角可行区间 =====\n");
    if (intervals.empty()) {
        printf("无可行臂角区间\n");
        return 0;
    }
    for (auto& iv : intervals)
        printf("  [%.4f, %.4f] rad  =  [%.2f, %.2f] deg\n",
               iv.first, iv.second, iv.first*180/M_PI, iv.second*180/M_PI);

    // ===== 反向验证：在每个区间边界处计算关节角，确认哪个关节到达限位 =====
    printf("\n===== 边界验证 (哪个关节在 ±175° 处) =====\n");

    // 复用 compute_arm_angle_intervals 内部的数学，内联一个小helper
    auto eval_at_psi = [&](double psi) {
        Eigen::Matrix3d R_0_d; Eigen::Vector3d P_0_d;
        getPoseFromArray(T_target, R_0_d, P_0_d);
        Eigen::Vector3d l_bs(0,0,dh[0][2]);
        Eigen::Vector3d l_wt(0,0,sqrt(dh[6][2]*dh[6][2]+dh[5][2]*dh[5][2]));
        Eigen::Vector3d x_sw = P_0_d - l_bs - R_0_d*l_wt;
        double norm_sw = x_sw.norm();
        Eigen::Vector3d u_sw = Eigen::Vector3d::Zero();
        if (norm_sw > 1e-12) u_sw = x_sw / norm_sw;

        double cos4 = std::clamp((norm_sw*norm_sw - dh[2][2]*dh[2][2] - dh[4][2]*dh[4][2])
                                  / (2.0*dh[2][2]*dh[4][2]), -1.0, 1.0);
        double theta_4 = -acos(cos4);
        double c4=cos(theta_4), s4=sin(theta_4);
        Eigen::Matrix3d R_3_4; R_3_4 << c4,0,-s4, 0,1,0, s4,0,c4;

        Eigen::Vector3d p_ref = Eigen::Vector3d(0,0,dh[2][2]) + R_3_4*Eigen::Vector3d(0,0,dh[4][2]);
        double px=p_ref(0),py=p_ref(1),pz=p_ref(2);
        double x=x_sw(0),y=x_sw(1),z=x_sw(2);
        double r=sqrt(px*px+pz*pz), t1r=0, t2r=0;
        if (r > 1e-12) {
            double alpha=atan2(px,pz), beta=acos(std::clamp(z/r,-1.0,1.0));
            auto Rz=[](double t){Eigen::Matrix3d R; double c=cos(t),s=sin(t); R<<c,-s,0,s,c,0,0,0,1; return R;};
            auto Ry=[](double t){Eigen::Matrix3d R; double c=cos(t),s=sin(t); R<<c,0,-s,0,1,0,s,0,c; return R;};
            double t2a=alpha+beta, t1a=atan2(y,x)-atan2(py,cos(t2a)*px-sin(t2a)*pz);
            double t2b=alpha-beta, t1b=atan2(y,x)-atan2(py,cos(t2b)*px-sin(t2b)*pz);
            double ea=(Rz(t1a)*Ry(t2a)*p_ref-x_sw).norm();
            double eb=(Rz(t1b)*Ry(t2b)*p_ref-x_sw).norm();
            if(ea<=eb){t1r=t1a;t2r=t2a;}else{t1r=t1b;t2r=t2b;}
        }
        auto Rz=[](double t){Eigen::Matrix3d R; double c=cos(t),s=sin(t); R<<c,-s,0,s,c,0,0,0,1; return R;};
        auto Ry=[](double t){Eigen::Matrix3d R; double c=cos(t),s=sin(t); R<<c,0,-s,0,1,0,s,0,c; return R;};
        Eigen::Matrix3d R03ref = Rz(t1r)*Ry(t2r);
        Eigen::Matrix3d ux; ux<<0,-u_sw(2),u_sw(1),u_sw(2),0,-u_sw(0),-u_sw(1),u_sw(0),0;
        Eigen::Matrix3d As=ux*R03ref, Bs=-ux*ux*R03ref, Cs=(u_sw*u_sw.transpose())*R03ref;
        Eigen::Matrix3d R03 = As*sin(psi)+Bs*cos(psi)+Cs;

        double s2=sqrt(R03(2,0)*R03(2,0)+R03(2,1)*R03(2,1));
        double t1,t2,t3;
        if(s2>1e-9){
            t2=atan2(s2,R03(2,2)); t1=atan2(-R03(1,2),-R03(0,2)); t3=atan2(-R03(2,1),R03(2,0));
        } else {
            t2=0; t3=0; t1=atan2(R03(1,0),R03(0,0));
        }

        Eigen::Matrix3d R47 = R_3_4.transpose()*R03.transpose()*R_0_d;
        double s6=sqrt(R47(2,0)*R47(2,0)+R47(2,1)*R47(2,1));
        double t5,t6,t7;
        if(s6>1e-9){
            t6=atan2(s6,R47(2,2)); t5=atan2(-R47(1,2),-R47(0,2)); t7=atan2(-R47(2,1),R47(2,0));
        } else {
            t6=0; t7=0; t5=atan2(R47(1,0),R47(0,0));
        }

        double q[7] = {t1*180/M_PI, t2*180/M_PI, t3*180/M_PI,
                       -theta_4*180/M_PI, t5*180/M_PI, t6*180/M_PI, t7*180/M_PI};
        printf("  psi=%7.2f°: J0=%7.1f J1=%7.1f J2=%7.1f J3=%7.1f J4=%7.1f J5=%7.1f J6=%7.1f",
               psi*180/M_PI, q[0],q[1],q[2],q[3],q[4],q[5],q[6]);
        for(int i=0;i<7;i++) if(fabs(q[i])>173.0) printf(" *** J%d@LIMIT",i);
        printf("\n");
    };

    for (auto& iv : intervals) {
        eval_at_psi(iv.first);   // 下边界
        eval_at_psi(iv.second);  // 上边界
    }
    return 0;
}
