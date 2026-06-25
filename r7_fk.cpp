#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double deg2rad(double deg) { return deg * M_PI / 180.0; }

void print_matrix(double T[4][4], const char* name) {
    cout << name << ":" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << setw(10) << setprecision(4) << T[i][j] << " ";
        }
        cout << endl;
    }
    cout << "位置 (mm): (" << T[0][3] << ", " << T[1][3] << ", " << T[2][3] << ")" << endl;
    cout << endl;
}

void print_matrix_cpp_format(double T[4][4]) {
    cout << "{";
    for (int i = 0; i < 4; i++) {
        cout << "{";
        for (int j = 0; j < 4; j++) {
            cout << setprecision(10) << T[i][j];
            if (j < 3) cout << ", ";
        }
        cout << "}";
        if (i < 3) cout << ",\n ";
    }
    cout << "}" << endl;
}

void print_pose(double T[4][4]) {
    cout << "\n=== 末端位姿 ===" << endl;
    cout << "位置 (mm): (" << T[0][3] << ", " << T[1][3] << ", " << T[2][3] << ")" << endl;

    // 提取旋转矩阵
    double r11 = T[0][0], r12 = T[0][1], r13 = T[0][2];
    double r21 = T[1][0], r22 = T[1][1], r23 = T[1][2];
    double r31 = T[2][0], r32 = T[2][1], r33 = T[2][2];

    // 旋转矩阵 → 四元数 (Shepperd's method)
    double q_w = sqrt(fmax(0, 1 + r11 + r22 + r33)) / 2;
    double q_x = sqrt(fmax(0, 1 + r11 - r22 - r33)) / 2;
    double q_y = sqrt(fmax(0, 1 - r11 + r22 - r33)) / 2;
    double q_z = sqrt(fmax(0, 1 - r11 - r22 + r33)) / 2;
    q_x = copysign(q_x, r32 - r23);
    q_y = copysign(q_y, r13 - r31);
    q_z = copysign(q_z, r21 - r12);

    cout << "四元数 (x,y,z,w): (" << q_x << ", " << q_y << ", " << q_z << ", " << q_w << ")" << endl;

    // 旋转矩阵 → RPY (ZYX Euler)
    double beta = atan2(-r31, sqrt(r11*r11 + r21*r21));
    double alpha, gamma;

    if (fabs(cos(beta)) > 1e-6) {
        alpha = atan2(r32/cos(beta), r33/cos(beta));
        gamma = atan2(r21/cos(beta), r11/cos(beta));
    } else {
        alpha = 0;
        gamma = atan2(r12, r11);
    }

    cout << "RPY (度): Roll=" << alpha*180/M_PI
         << ", Pitch=" << beta*180/M_PI
         << ", Yaw=" << gamma*180/M_PI << endl;
}

void mdh_transform(double alpha, double a, double d, double theta, double T[4][4]) {
    double ca = cos(deg2rad(alpha));
    double sa = sin(deg2rad(alpha));
    double ct = cos(deg2rad(theta));
    double st = sin(deg2rad(theta));
    
    T[0][0] = ct;      T[0][1] = -st;     T[0][2] = 0;   T[0][3] = a;
    T[1][0] = st*ca;   T[1][1] = ct*ca;   T[1][2] = -sa; T[1][3] = -sa*d;
    T[2][0] = st*sa;   T[2][1] = ct*sa;   T[2][2] = ca;  T[2][3] = ca*d;
    T[3][0] = 0;       T[3][1] = 0;       T[3][2] = 0;   T[3][3] = 1;
}

void multiply(double A[4][4], double B[4][4], double result[4][4]) {
    double temp[4][4] = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                temp[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            result[i][j] = temp[i][j];
}

void forward_kinematics(double joint_angles[7], bool verbose = false) {
    // DH参数表: [αᵢ₋₁, aᵢ₋₁, dᵢ, θᵢ_offset]
    // 实际θᵢ = θᵢ_offset + joint_angles[i]
    double dh[7][4] = {
        {   0,  0,   84+95,  0 },
        {  -90,  0,   0,   0 },
        { 90,  0,  215+260,  0 },  // 84+215
        {  -90,  0,   0,   0 },
        { 90,  0,  415+60,  0 },  // 260+415
        {  -90,  0, -0.0,  0 },
        { 90,  0,  145,  0 },  // 60+145
    };
    
    double T_total[4][4] = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}};
    double Ti[4][4], T_temp[4][4];
    
    if (verbose) {
        cout << "=== 输入关节角度 (度) ===" << endl;
        for (int i = 0; i < 7; i++) {
            cout << "关节" << i+1 << ": " << joint_angles[i] << "°" << endl;
        }
        cout << "\n=== 逐轴累积变换 ===" << endl;
    }
    
    for (int i = 0; i < 7; i++) {
        // 实际角度 = DH表中的偏置 + 输入角度
        double theta = dh[i][3] + joint_angles[i];
        mdh_transform(dh[i][0], dh[i][1], dh[i][2], theta, Ti);

        multiply(T_total, Ti, T_temp);
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
                T_total[r][c] = T_temp[r][c];

        if (verbose) {
            cout << "T0_" << i+1 << " = double [4][4] = ";
            print_matrix_cpp_format(T_total);
            cout << endl;
        }
    }
    
    print_pose(T_total);
}

int main() {
    double joint_angles[7];
    
    cout << "=== 7自由度机械臂正运动学计算 ===" << endl;
    cout << "请输入7个关节角度 (度)，用空格分隔:" << endl;
    cout << "格式: θ1 θ2 θ3 θ4 θ5 θ6 θ7" << endl;
    cout << "(输入 0 0 0 0 0 0 0 可验证零位)" << endl;
    cout << "> ";
    
    for (int i = 0; i < 7; i++) {
        cin >> joint_angles[i];
    }
    
    // 询问是否显示详细过程
    char show_detail;
    cout << "\n是否显示详细变换过程? (y/n): ";
    cin >> show_detail;
    
    bool verbose = (show_detail == 'y' || show_detail == 'Y');
    
    forward_kinematics(joint_angles, verbose);
    
    // 验证零位
    if (joint_angles[0] == 0 && joint_angles[1] == 0 && joint_angles[2] == 0 &&
        joint_angles[3] == 0 && joint_angles[4] == 0 && joint_angles[5] == 0 &&
        joint_angles[6] == 0) {
        cout << "\n=== 零位验证 ===" << endl;
        cout << "预期位置: (0, -2.5, 1274) mm" << endl;
        double z_sum = 95 + 84 + 215 + 260 + 415 + 60 + 145;
        cout << "理论 z 高度: " << z_sum << " mm" << endl;
    }
    
    return 0;
}
