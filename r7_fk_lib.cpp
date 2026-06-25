#include <cmath>

static double deg2rad(double deg) { return deg * M_PI / 180.0; }

static void mdh_transform(double alpha, double a, double d, double theta, double T[4][4]) {
    double ca = cos(deg2rad(alpha));
    double sa = sin(deg2rad(alpha));
    double ct = cos(deg2rad(theta));
    double st = sin(deg2rad(theta));

    T[0][0] = ct;      T[0][1] = -st;     T[0][2] = 0;   T[0][3] = a;
    T[1][0] = st*ca;   T[1][1] = ct*ca;   T[1][2] = -sa; T[1][3] = -sa*d;
    T[2][0] = st*sa;   T[2][1] = ct*sa;   T[2][2] = ca;  T[2][3] = ca*d;
    T[3][0] = 0;       T[3][1] = 0;       T[3][2] = 0;   T[3][3] = 1;
}

static void multiply(double A[4][4], double B[4][4], double result[4][4]) {
    double temp[4][4] = {0};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                temp[i][j] += A[i][k] * B[k][j];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            result[i][j] = temp[i][j];
}

// ============================================================
// 正运动学
//   joint_angles[7] — 关节角度 (度)
//   T_out[4][4]     — 输出齐次矩阵 (位置单位: mm), 无需预初始化
//   joint           — 输出第几个关节的位姿 (1~7, 默认 7 = 末端)
// ============================================================
void forward_kinematics(const double joint_angles[7], double T_out[4][4], int joint) {
    double dh[7][4] = {
        {   0,  0,   84+95,   0 },
        {  -90,  0,   0,       0 },
        {  90,  0,  215+260,   0 },
        {  -90,  0,   0,       0 },
        {  90,  0,  415+60,    0 },
        {  -90,  0, -2.5,       0 },
        {  90,  0,  145,        0 },
    };

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            T_out[i][j] = (i == j) ? 1.0 : 0.0;

    double Ti[4][4], T_temp[4][4];
    for (int i = 0; i < joint; i++) {
        double theta = dh[i][3] + joint_angles[i];
        mdh_transform(dh[i][0], dh[i][1], dh[i][2], theta, Ti);
        multiply(T_out, Ti, T_temp);
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
                T_out[r][c] = T_temp[r][c];
    }
}
