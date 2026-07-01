#include <iostream>
#include <cmath>
#include <iomanip>
#include <Eigen/Geometry>

void forward_kinematics(const double joint_angles[7], double T_out[4][4], int joint = 7);

const double dh[7][4] = {
    {   0,  0,   84+95,   0 },
    {  -90,  0,   0,       0 },
    {  90,  0,  215+260,   0 },
    {  -90,  0,   0,       0 },
    {  90,  0,  415+60,    0 },
    {  -90,  0, -2.5,       0 },
    {  90,  0,  145,        0 },
};

void getPoseFromArray(const double T[4][4], Eigen::Matrix3d& R, Eigen::Vector3d& p) {
    R << T[0][0], T[0][1], T[0][2],
         T[1][0], T[1][1], T[1][2],
         T[2][0], T[2][1], T[2][2];
    p << T[0][3], T[1][3], T[2][3];
}

int main() {
    std::cout << std::fixed << std::setprecision(4);

    // 测试几个姿态
    double test_cases[3][7] = {
        {0, 0, 0, 0, 0, 0, 0},           // 零位
        {30, 45, 60, 30, 45, 60, 30},    // 随机姿态1
        {-30, 60, -45, 90, 30, -60, 45}  // 随机姿态2
    };

    for (int test = 0; test < 3; test++) {
        std::cout << "\n========== 测试姿态 " << test+1 << " ==========\n";
        std::cout << "关节角度: [";
        for (int i = 0; i < 7; i++) std::cout << test_cases[test][i] << (i<6 ? ", " : "]\n");

        // 计算 FK
        double T_7[4][4], T_6[4][4];
        forward_kinematics(test_cases[test], T_7, 7);  // 末端
        forward_kinematics(test_cases[test], T_6, 6);  // 腕中心

        Eigen::Matrix3d R_7; Eigen::Vector3d P_7;
        getPoseFromArray(T_7, R_7, P_7);

        Eigen::Vector3d P_6_actual(T_6[0][3], T_6[1][3], T_6[2][3]);

        // 旧方法（错误）：sqrt(145^2 + 2.5^2) = 145.02
        Eigen::Vector3d l_wt_old(0, 0, sqrt(dh[6][2]*dh[6][2] + dh[5][2]*dh[5][2]));
        Eigen::Vector3d P_6_old = P_7 - R_7 * l_wt_old;

        // 新方法（正确）：只用 -145
        Eigen::Vector3d l_wt_new(0, 0, -dh[6][2]);
        Eigen::Vector3d P_6_new = P_7 - R_7 * l_wt_new;

        std::cout << "\n末端位置 P_7:          (" << P_7.transpose() << ")\n";
        std::cout << "实际腕中心 P_6 (FK):   (" << P_6_actual.transpose() << ")\n";
        std::cout << "旧方法计算的 P_6:       (" << P_6_old.transpose() << ")\n";
        std::cout << "新方法计算的 P_6:       (" << P_6_new.transpose() << ")\n";

        double error_old = (P_6_old - P_6_actual).norm();
        double error_new = (P_6_new - P_6_actual).norm();

        std::cout << "\n旧方法误差: " << error_old << " mm\n";
        std::cout << "新方法误差: " << error_new << " mm "
                  << (error_new < 1e-6 ? "✓ 正确" : "✗ 错误") << "\n";
    }

    return 0;
}
