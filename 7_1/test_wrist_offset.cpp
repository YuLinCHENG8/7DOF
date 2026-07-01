#include <iostream>
#include <cmath>
#include <iomanip>

void forward_kinematics(const double joint_angles[7], double T_out[4][4], int joint = 7);

int main() {
    // 测试零位姿态
    double q_zero[7] = {0, 0, 0, 0, 0, 0, 0};

    double T_5[4][4], T_6[4][4], T_7[4][4];
    forward_kinematics(q_zero, T_5, 5);  // Frame 5 (关节5后)
    forward_kinematics(q_zero, T_6, 6);  // Frame 6 (关节6后，腕中心)
    forward_kinematics(q_zero, T_7, 7);  // Frame 7 (末端)

    std::cout << std::fixed << std::setprecision(2);

    std::cout << "Frame 5 位置: (" << T_5[0][3] << ", " << T_5[1][3] << ", " << T_5[2][3] << ")\n";
    std::cout << "Frame 6 位置: (" << T_6[0][3] << ", " << T_6[1][3] << ", " << T_6[2][3] << ")\n";
    std::cout << "Frame 7 位置: (" << T_7[0][3] << ", " << T_7[1][3] << ", " << T_7[2][3] << ")\n";

    // 计算从 Frame 7 到 Frame 6 的向量（在 Frame 7 坐标系下）
    double dx = T_6[0][3] - T_7[0][3];
    double dy = T_6[1][3] - T_7[1][3];
    double dz = T_6[2][3] - T_7[2][3];

    std::cout << "\n从末端到腕中心的向量（在 base frame 下）:\n";
    std::cout << "  Δx = " << dx << " mm\n";
    std::cout << "  Δy = " << dy << " mm\n";
    std::cout << "  Δz = " << dz << " mm\n";
    std::cout << "  长度 = " << std::sqrt(dx*dx + dy*dy + dz*dz) << " mm\n";

    // 提取 Frame 7 的旋转矩阵
    std::cout << "\nFrame 7 旋转矩阵:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "  [";
        for (int j = 0; j < 3; j++) {
            std::cout << std::setw(8) << T_7[i][j];
        }
        std::cout << " ]\n";
    }

    // 在 Frame 7 坐标系下，从末端到腕中心的向量
    // T_6 = T_7 * T_7_6^(-1)
    // 所以 p_6 = p_7 + R_7 * p_7_6
    // 反过来：p_7_6 = R_7^T * (p_6 - p_7)

    double R7_T[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R7_T[i][j] = T_7[j][i];  // 转置

    double p_local[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            p_local[i] += R7_T[i][j] * (j == 0 ? dx : j == 1 ? dy : dz);

    std::cout << "\n从末端到腕中心的向量（在 Frame 7 局部坐标系下）:\n";
    std::cout << "  x_local = " << p_local[0] << " mm\n";
    std::cout << "  y_local = " << p_local[1] << " mm\n";
    std::cout << "  z_local = " << p_local[2] << " mm\n";

    return 0;
}
