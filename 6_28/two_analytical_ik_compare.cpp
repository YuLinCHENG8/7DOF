#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <utility> 
#include <cmath>
#include <chrono>
#include <Eigen/Geometry>
#include <iomanip>
#include <cstring>

using std::vector;
using std::pair;
using namespace std;

#define CONTROL_HZ  500.0
#define DT          (1.0 / CONTROL_HZ)
#define N_JOINTS    7
#define N_PSI       7

//璁＄畻鍑虹殑鑷傝鍙鍖洪棿
constexpr double ARM_ANGLE_DISTRICT_1_1 = 0.0;
constexpr double ARM_ANGLE_DISTRICT_1_2 = 0.9021;
constexpr double ARM_ANGLE_DISTRICT_2_1 = 1.1019;
constexpr double ARM_ANGLE_DISTRICT_2_2 = 2.646297;

// 鍏宠妭闄愪綅锛堝姬搴︼級175*M_PI/180 = 3.054326
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


// 璁＄畻涓や釜榻愭鍙樻崲鐭╅樀鐨勫樊寮?
void compareTransforms(const double T1[4][4], const double T2[4][4]) {
    // 鎻愬彇骞崇Щ閮ㄥ垎 (鍋囪鐭╅樀鎸夎涓诲簭瀛樺偍锛孴[row][col])
    double tx1 = T1[0][3], ty1 = T1[1][3], tz1 = T1[2][3];
    double tx2 = T2[0][3], ty2 = T2[1][3], tz2 = T2[2][3];
    
    // 浣嶇疆璇樊锛氬钩绉诲悜閲忕殑娆у嚑閲屽緱璺濈
    double dx = tx1 - tx2;
    double dy = ty1 - ty2;
    double dz = tz1 - tz2;
    double position_error = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // 鎻愬彇鏃嬭浆鐭╅樀閮ㄥ垎 (3x3)
    double R1[3][3], R2[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            R1[i][j] = T1[i][j];
            R2[i][j] = T2[i][j];
        }
    }
    
    // 璁＄畻鐩稿鏃嬭浆鐭╅樀 螖R = R1 * R2^T
    // 棣栧厛璁＄畻 R2 鐨勮浆缃?R2_T
    double R2_T[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            R2_T[i][j] = R2[j][i];
    
    // 鐭╅樀涔樻硶 螖R = R1 * R2_T
    double delta_R[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                delta_R[i][j] += R1[i][k] * R2_T[k][j];
    
    // 閫氳繃杩?鍙嶅绉伴儴鍒嗚绠楁棆杞搴︼紝鐢?atan2 瀵瑰皬瑙掑害鏇寸ǔ瀹?
    double trace = delta_R[0][0] + delta_R[1][1] + delta_R[2][2];
    double s = std::sqrt((delta_R[2][1] - delta_R[1][2]) * (delta_R[2][1] - delta_R[1][2]) +
                         (delta_R[0][2] - delta_R[2][0]) * (delta_R[0][2] - delta_R[2][0]) +
                         (delta_R[1][0] - delta_R[0][1]) * (delta_R[1][0] - delta_R[0][1])) / 2.0;
    double c = (trace - 1.0) / 2.0;
    double angle_error_rad = std::atan2(s, c);  // 寮у害
    double angle_error_deg = angle_error_rad * 180.0 / M_PI;  // 杞崲涓哄害
    
    // 杈撳嚭瀵规瘮缁撴灉
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "(Position error): " << position_error << " units\n";
    std::cout << "(Orientation error): " << angle_error_rad << " rad ("
              << angle_error_deg << " deg)\n";
}


//鏂板缓鍙嶅鎴愮煩闃?
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


double arm_plane_angle_test(const double q[7]) {
    Eigen::Matrix3d R,R_0_4,R_0_4_ref;
    Eigen::Vector3d E, E_ref;

    fK_eigen(q, R_0_4, E, 4);  // elbow 

    double q0[7];
    for (int i = 0; i < 7; i++) q0[i] = q[i];
    q0[2] = 0.0;           // theta_3 = 0
    fK_eigen(q0, R_0_4_ref, E_ref, 4); // elbow ref (涓嶉殢 theta_3 鍙樺寲)

    //--------------鍏堣绠?R_0_arm_angle
    Eigen::Matrix3d R_0_arm_angle;
    R_0_arm_angle = R_0_4 * R_0_4_ref.transpose();

    //------------璁＄畻 x_0_sw锛屽€熼壌涓嬮潰ik鐨勮绠楅€昏緫
    Eigen::Matrix3d R_0_desire;
    Eigen::Vector3d P_0_desire;
    //鑾峰彇base涓嬬殑SW鍧愭爣,鍙婂叾鍗曚綅鍚戦噺
    Eigen::Vector3d x_0_sw; // 鍦?base frame涓?鑾峰彇sw鍚戦噺
    Eigen::Vector3d u_0_sw;//x_sw_0 鐨勫崟浣嶅悜閲?
    double T_target[4][4] = {};
    forward_kinematics(q,T_target,7);

    getPoseFromArray(T_target,R_0_desire,P_0_desire);//鑾峰彇鐩爣浣嶅Э鐨刾osi rot
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
    double cos_psi = 0.5 * (R_0_arm_angle.trace() - 1.0);

    if (cos_psi > 1.0) cos_psi = 1.0;
    if (cos_psi < -1.0) cos_psi = -1.0;

    Eigen::Vector3d vee;
    vee << R_0_arm_angle(2, 1) - R_0_arm_angle(1, 2),
           R_0_arm_angle(0, 2) - R_0_arm_angle(2, 0),
           R_0_arm_angle(1, 0) - R_0_arm_angle(0, 1);

    double sin_psi = 0.5 * u_0_sw.dot(vee);

    double arm_angle = std::atan2(sin_psi, cos_psi);

    printf("arm_angle = %f rad, %f deg\n",
           arm_angle,
           arm_angle * 180.0 / M_PI);

    return arm_angle;
}

void analytical_ik_test(const double T_target[4][4]/* , const double q_init[7], double psi, double q_out[7] */){

    double arm_angle = 0;

    Eigen::Matrix3d R_0_desire;
    Eigen::Vector3d P_0_desire;
    //鑾峰彇base涓嬬殑SW鍧愭爣,鍙婂叾鍗曚綅鍚戦噺
    Eigen::Vector3d x_0_sw; // 鍦?base frame涓?鑾峰彇sw鍚戦噺
    Eigen::Vector3d u_0_sw;//x_sw_0 鐨勫崟浣嶅悜閲?

    getPoseFromArray(T_target,R_0_desire,P_0_desire);//鑾峰彇鐩爣浣嶅Э鐨刾osi rot
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

    //==============瑙?theta 4=========== checked
    double cos_theta4 = (pow(norm_x_0_sw,2) - pow(dh[2][2],2) - pow(dh[4][2],2)) / ( 2* dh[2][2] * dh[4][2]);
    
    if (cos_theta4 < -1.0) cos_theta4 = -1.0;
    if (cos_theta4 > 1.0) cos_theta4 = 1.0;
    double theta_4 = -std::acos(cos_theta4);
    printf("------------第四个关节角度为：%f \n",theta_4*180/M_PI);

    //閫氳繃缃楀痉閲屾牸鏂彉鎹㈡眰鍑篟_0_3锛岃繘鑰屾眰鍑?theta0 - theta3
    Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d ux = skew(u_0_sw);
    Eigen::Matrix3d R_0_armangle =I3 + std::sin(arm_angle) * ux + (1.0 - std::cos(arm_angle)) * ux * ux;
    double theta_1[2] = {};//+ -
    double theta_2[2] = {};// + -
    double theta_3[2] = {};// + -
    double theta_5[2] = {};
    double theta_6[2] = {};
    double theta_7[2] = {};
    double cos_theta2 = - R_0_armangle(2,1);
    theta_2[0] = acos(cos_theta2);
    theta_2[1] = -acos(cos_theta2);
    Eigen::Matrix3d R_4_7[2];
    Eigen::Matrix3d R_0_4[2];
    Eigen::Vector3d posi_temp[2];
    for(int i = 0;i<2;i++){
        double s2 = sin(theta_2[i]);
        if (std::abs(s2) < 1e-6) {
            theta_1[i] = atan2(-R_0_armangle(1,1), -R_0_armangle(0,1));
            theta_3[i] = 0.0;
        } else {
            theta_1[i] = atan2(-R_0_armangle(1,1)/s2, -R_0_armangle(0,1)/s2);
            theta_3[i] = atan2(-R_0_armangle(2,2)/s2, -R_0_armangle(2,0)/s2);
        }

        double joint_angle[7] = {theta_1[i],theta_2[i],theta_3[i],theta_4,0,0,0};
        double T_0_4[4][4] = {};
        forward_kinematics(joint_angle,T_0_4,4);
        getPoseFromArray(T_0_4,R_0_4[i],posi_temp[i]);
        R_4_7[i] = R_0_4[i].transpose() * R_0_desire;
        double r13 = R_4_7[i](0, 2);
        double r23 = R_4_7[i](1, 2);
        double r31 = R_4_7[i](2, 0);
        double r32 = R_4_7[i](2, 1);
        double r33 = R_4_7[i](2, 2);



        theta_6[0] = acos(-r33);
        theta_6[1] = -acos(-r33);

        theta_5[i] = std::atan2(r23/sin(theta_6[i]), r13/sin(theta_6[i]));

        theta_7[i] = std::atan2(r32/sin(theta_6[i]), -r31/sin(theta_6[i]));
        printf("第 %d 组解的theta_1为:%f, theta_2为:%f, theta_3为:%f,theta_4为:%f, theta_5为:%f, theta_6为:%f , theta_7为:%f\n",
            i+1 ,theta_1[i]*180/M_PI, theta_2[i]*180/M_PI, theta_3[i]*180/M_PI,theta_4*180/M_PI,
            theta_5[i]*180/M_PI, theta_6[i]*180/M_PI, theta_7[i]*180/M_PI);


    }
    


    

    

}



/**
 * @brief 璁＄畻涓や釜骞抽潰鐨勫す瑙掞細
 *       鍙傝€冨钩闈細theta_3 = 0 鏃讹紝shoulder(frame2)銆乪lbow(frame3)銆亀rist(frame5) 纭畾鐨勫钩闈?
 *       瀹為檯骞抽潰锛歵heta_3 鈮?0 鏃剁殑鍚屼竴骞抽潰
 * @param q 褰撳墠瑙掑害
 * @return psi 杩斿洖寮у害
 */
double arm_plane_angle(const double q[7]) {
    // 1. 鑾峰彇涓変釜鍏抽敭鐐笶igen::Matrix3d R_elbow; Eigen::Vector3d p_elbow;
    Eigen::Matrix3d R_elbow; 
    Eigen::Vector3d p_elbow;
    fK_eigen(q, R_elbow, p_elbow, 4);  // 鑲樺叧鑺備綅缃?

    double T[4][4] = {};
    forward_kinematics(q, T, 7);
    Eigen::Matrix3d R_ee; 
    Eigen::Vector3d p_ee;
    getPoseFromArray(T, R_ee, p_ee);

    // 2. 鑲╀腑蹇冿紙鍥哄畾锛?
    Eigen::Vector3d p_s(0, 0, dh[0][2]);

    // 3. 鑵曚腑蹇?
    Eigen::Vector3d l_7_wt(0, 0, sqrt(pow(dh[6][2],2) + pow(dh[5][2],2)));
    Eigen::Vector3d p_w = p_ee - R_ee * l_7_wt;

    // 4. 鑲┾啋鑵曞崟浣嶅悜閲忥紙鏃嬭浆杞达級
    Eigen::Vector3d sw = p_w - p_s;
    Eigen::Vector3d u = sw.normalized();

    // 5. 鍙傝€冩柟鍚戯細涓栫晫Z杞存姇褰卞埌鍨傜洿u鐨勫钩闈?
    Eigen::Vector3d z_world(0, 0, 1);
    Eigen::Vector3d v_ref = z_world - z_world.dot(u) * u;
    if (v_ref.norm() < 1e-6) {
        // sw杩戜技骞宠Z杞存椂鏀圭敤X杞?
        Eigen::Vector3d x_world(1, 0, 0);
        v_ref = x_world - x_world.dot(u) * x_world;
    }
    v_ref.normalize();

    // 6. 褰撳墠鑲樺悜閲忔姇褰卞埌鍨傜洿u鐨勫钩闈?
    Eigen::Vector3d se = p_elbow - p_s;
    Eigen::Vector3d v_e = se - se.dot(u) * u;
    if (v_e.norm() < 1e-6) return 0.0;  // 濂囧紓锛氳倶鍦ㄨ偐鑵曡繛绾夸笂
    v_e.normalize();

    // 7. 甯︾鍙疯搴?
    double cos_psi = std::clamp(v_ref.dot(v_e), -1.0, 1.0);
    double sin_psi = u.dot(v_ref.cross(v_e));
    return std::atan2(sin_psi, cos_psi);
}


/**
 * @brief 鍗曡皟鍨?(Delta3 < 0) 姹傝В鑷傝鍙鍖洪棿
 * 
 * 鍗曡皟鍨嬬壒鐐癸細
 *   - 胃(蠄) 鍦ㄦ暣涓?[-蟺, 蟺) 涓婂崟璋冮€掑鎴栭€掑噺
 *   - 鏃犳瀬鍊肩偣
 *   - 閫氳繃 胃_min/胃_max 鐩存帴瑙ｅ搴旂殑 蠄锛屽啀鏍规嵁鍗曡皟鏂瑰悜纭畾鍖洪棿
 * 
 * @param A_n, B_n, C_n    鍒嗗瓙绯绘暟: num = A_n*sin蠄 + B_n*cos蠄 + C_n
 * @param A_d, B_d, C_d    鍒嗘瘝绯绘暟: den = A_d*sin蠄 + B_d*cos蠄 + C_d
 * @param a_t, b_t, c_t    瀵兼暟绯绘暟锛堢敤浜庣‘瀹氬崟璋冩柟鍚戯級
 * @param theta_min        鍏宠妭瑙掍笅闄愶紙寮у害锛?
 * @param theta_max        鍏宠妭瑙掍笂闄愶紙寮у害锛?
 * @param intervals        杈撳嚭锛氬彲琛岃噦瑙掑尯闂村垪琛?
 * @return true  姹傝В鎴愬姛
 * @return false 姹傝В澶辫触锛堝鏃犺В锛?
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
    // 1. 纭畾鍗曡皟鏂瑰悜
    //    鍦?蠄=0 澶勶紝瀵兼暟鐨勫垎瀛愪负: a_t*sin(0) + b_t*cos(0) + c_t = b_t + c_t
    //    鍒嗘瘝鎭掓锛屾墍浠ュ鏁扮殑绗﹀彿鐢?(b_t + c_t) 鍐冲畾
    // ============================================================
    double sign_val = b_t + c_t;
    bool increasing = (sign_val > 0.0);   // true: 胃闅徬堥€掑; false: 閫掑噺

    // ============================================================
    // 2. 杈呭姪鍑芥暟锛氬浜庣粰瀹氱殑 tan_target锛岃В瀵瑰簲鐨?蠄锛堝綊涓€鍖栧埌 [-蟺, 蟺)锛?
    //    鏂圭▼锛?(A_n sin蠄 + B_n cos蠄 + C_n) / (A_d sin蠄 + B_d cos蠄 + C_d) = tan_target
    // ============================================================
    auto solve_psi_for_tan = [&](double tan_target) -> double {
        double A = A_n - tan_target * A_d;
        double B = B_n - tan_target * B_d;
        double C = C_n - tan_target * C_d;

        double r = sqrt(A * A + B * B);
        if (std::abs(C) > r) return NAN;   // 鏃犺В

        // 瑙?A sin蠄 + B cos蠄 = -C
        double phi = atan2(B, A);           // A sin蠄 + B cos蠄 = r sin(蠄 + phi)
        double sin_val = -C / r;
        double alpha = asin(sin_val);       // 蠄 + phi = alpha 鎴?蟺 - alpha

        double psi_1 = alpha - phi;
        double psi_2 = M_PI - alpha - phi;

        // 褰掍竴鍖栧埌 [-蟺, 蟺)
        auto wrap = [](double x) {
            while (x < -M_PI) x += 2 * M_PI;
            while (x >=  M_PI) x -= 2 * M_PI;
            return x;
        };
        psi_1 = wrap(psi_1);
        psi_2 = wrap(psi_2);

        // 妫€楠屽摢涓?蠄 鐪熺殑婊¤冻 tan(胃(蠄)) = tan_target
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
    // 3. 瀵?theta_min 鍜?theta_max 姹傚搴旂殑鑷傝杈圭晫
    // ============================================================
    double psi_L = solve_psi_for_tan(tan(theta_min));
    double psi_U = solve_psi_for_tan(tan(theta_max));

    if (std::isnan(psi_L) || std::isnan(psi_U)) {
        // 鏃犺В锛氶檺浣嶅唴娌℃湁鍙鑷傝
        return false;
    }

    // ============================================================
    // 4. 鏍规嵁鍗曡皟鏂瑰悜纭畾鍙鍖洪棿
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
    // 5. 杈撳嚭鍖洪棿锛堣€冭檻璺ㄨ秺 -蟺 鐨勬儏鍐碉級
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
 * @brief 闇€瑕佷紶鍏ヨ偐鍏宠妭鍜岃厱鍏宠妭涔嬮棿鐨勮繛绾?
 * @param q_init 涓鸿搴?
 * @return sigualrity_type 濂囧紓绫诲瀷
 */
sigualrity_type check_near_singularity(double q_init[7], Eigen::Vector3d x_0_sw){
    //1銆佽偐閮ㄥ寮傦紝鏉′欢涓簒_0_sw 脳 Z0 = 0銆?涓€鑸笉浼氬埌杈?
    
    // 鍩哄骇鏍囩郴 Z 杞?(0,0,1)
    Eigen::Vector3d z0(0.0, 0.0, 1.0);
    // 璁＄畻鍙夌Н
    Eigen::Vector3d cross_result = x_0_sw.cross(z0);
    if(fabs(cross_result.norm()) < 0.05){ //褰撹厱閮ㄨ窛绂粃杞寸嚎 5cm
        return shoulder_singular;
    }
    //2銆佽倶閮ㄥ寮傦紝璁＄畻鍑烘潵鐨?theta4 涓嶄负 0/pi銆?涓€鑸彧鏈夊伐浣滅┖闂磋竟鐣屾墠鏈?
    else if(fabs(q_init[3]) < 10.0 || fabs(q_init[3]+180) < 10.0 ||fabs(q_init[3]-180) < 10.0 ){
        return elbow_singular;
    }

    //3銆佽厱閮ㄥ寮?theta6 瑙掑害鎺ヨ繎 0锛屽叧鑺?5銆?6 鍏辩嚎鐨勬儏鍐?
    else if(fabs(q_init[5]) < 10.0 || fabs(q_init[5]+180) < 10.0 ||fabs(q_init[5]-180) < 10.0 ){
        return wrist_singular;
    }

    //4銆佹枃绔犱腑鎻愬強鐨?sin theta2锛? 涓洪浂鐨勬儏鍐?
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




// ============================================================
// 以下为从 tsinghua_paper.py 移植的几何求解版本 (函数以 _py 结尾)
// ============================================================

// 余弦定理求 theta4 (py: _solve_theta4_from_triangle)
//
static double solve_theta4_from_triangle_py(const Eigen::Vector3d& S, const Eigen::Vector3d& W)
{
    double l_sw = (W - S).norm();
    double l_se = std::abs(dh[2][2]);
    double l_ew = std::abs(dh[4][2]);
    double c4 = (l_sw*l_sw - l_se*l_se - l_ew*l_ew) / (2.0 * l_se * l_ew);
    if (c4 < -1.0) c4 = -1.0;
    if (c4 > 1.0) c4 = 1.0;
    return std::acos(c4);
}

// 从目标位姿求 S, W, theta4_abs, u_sw (py: _compute_swe_from_target)
static void compute_swe_from_target_py(
    const Eigen::Matrix3d& R_des, const Eigen::Vector3d& p_des,
    Eigen::Vector3d& S, Eigen::Vector3d& W,
    double& q4_abs, Eigen::Vector3d& u_sw)
{
    double d1_val = dh[0][2];
    double d7_val = dh[6][2];
    double d8_val = 0.0;

    Eigen::Vector3d z7 = R_des.col(2);

    // py: O7 = p_target - d8*z7 = p_target
    // py: W = O7 - d6*z7, 其中 d6 = p.d_i[6] = d7_in_dh = 145
    // py 漏掉了 -2.5 偏移；保持 py 原逻辑以对比
    W = p_des - d7_val * z7;

    S << 0.0, 0.0, d1_val;

    q4_abs = solve_theta4_from_triangle_py(S, W);
    printf("============py 版本的臂角数据为：%f \n",q4_abs*180/M_PI);

    Eigen::Vector3d v_sw = W - S;
    double n_sw = v_sw.norm();
    if (n_sw > 1e-12)
        u_sw = v_sw / n_sw;
    else
        u_sw = Eigen::Vector3d(0.0, 0.0, 1.0);
}

// 由臂角 theta0 求肘部位置 E (py: _elbow_from_arm_angle)
static Eigen::Vector3d elbow_from_arm_angle_py(
    const Eigen::Vector3d& S, const Eigen::Vector3d& W,
    double theta0)
{
    double l_se = std::abs(dh[2][2]);
    double l_ew = std::abs(dh[4][2]);
    Eigen::Vector3d sw = W - S;
    double l_sw = sw.norm();
    Eigen::Vector3d u_sw = sw / l_sw;

    double x = (l_se*l_se - l_ew*l_ew + l_sw*l_sw) / (2.0 * l_sw);
    double r2 = l_se*l_se - x*x;
    double r = std::sqrt(std::max(0.0, r2));
    Eigen::Vector3d C = S + x * u_sw;

    Eigen::Vector3d os_vec = S;
    Eigen::Vector3d t = os_vec.cross(u_sw);
    if (t.norm() < 1e-12) {
        t = Eigen::Vector3d(1.0, 0.0, 0.0);
        if (std::abs(t.dot(u_sw)) > 0.999) t = Eigen::Vector3d(0.0, 1.0, 0.0);
        t = t - t.dot(u_sw)*u_sw;
    }
    Eigen::Vector3d e1 = t.normalized();
    Eigen::Vector3d e2 = u_sw.cross(e1).normalized();

    Eigen::Vector3d E = C + r * (std::cos(theta0) * e1 + std::sin(theta0) * e2);
    return E;
}

// 从肘位置 E, 腕心 W, q4 计算 q1,q2,q3 (py: _solve_q123_from_swe)
static int solve_q123_from_swe_py(
    const Eigen::Vector3d& E, const Eigen::Vector3d& W, double q4,
    double sols[2][3])
{
    double d0 = dh[0][2];
    double d2 = dh[2][2];
    double d4_val = dh[4][2];
    double Ex=E.x(), Ey=E.y(), Ez=E.z();

    double c2 = (Ez - d0) / d2;
    if (c2 < -1.0) c2 = -1.0;
    if (c2 > 1.0) c2 = 1.0;
    double s2_abs = std::sqrt(std::max(0.0, 1.0 - c2*c2));

    double s4 = std::sin(q4);
    double c4 = std::cos(q4);

    int nsol = 0;
    for (int sign = 0; sign < 2; sign++) {
        double s2 = (sign == 0) ? s2_abs : -s2_abs;

        double c1 = -Ex / (d2 * s2);
        double s1 = -Ey / (d2 * s2);
        double n1 = std::hypot(c1, s1);
        if (n1 < 1e-15) continue;
        c1 /= n1; s1 /= n1;
        double q1 = std::atan2(s1, c1);
        double q2 = std::atan2(s2, c2);

        Eigen::Vector3d v = W - E;
        Eigen::Vector3d col2 = -v / d4_val;
        double u1 = col2.x(), u2 = col2.y(), u3 = col2.z();

        double b1 = (s2*c1*c4 - u1) / s4;
        double b2 = (u2 - s1*s2*c4) / s4;
        double s3 = s1*b1 + c1*b2;
        double c2c3 = -c1*b1 + s1*b2;
        double c3;
        if (std::abs(c2) > 1e-8) {
            c3 = c2c3 / c2;
        } else {
            c3 = (u3 + c2*c4) / (s2*s4);
        }
        double n3 = std::hypot(s3, c3);
        if (n3 < 1e-15) continue;
        s3 /= n3; c3 /= n3;
        double q3 = std::atan2(s3, c3);

        sols[nsol][0] = q1;
        sols[nsol][1] = q2;
        sols[nsol][2] = q3;
        nsol++;
    }
    return nsol;
}

// 从 T47 提取 wrist: q5,q6,q7 (py: _extract_567_from_T47_paper)
static int extract_567_from_T47_py(const Eigen::Matrix3d& T47, double sols[2][3])
{
    int nsol = 0;
    double c6 = T47(1, 2);
    if (c6 < -1.0) c6 = -1.0;
    if (c6 > 1.0) c6 = 1.0;
    for (int sign = 0; sign < 2; sign++) {
        double s6 = (sign == 0 ? 1.0 : -1.0) * std::sqrt(std::max(0.0, 1.0 - c6*c6));
        if (std::abs(s6) < 1e-8) continue;
        double th6 = std::atan2(s6, c6);
        double th5 = std::atan2(T47(2, 2) / s6, T47(0, 2) / s6);
        double th7 = std::atan2(T47(1, 1) / s6, -T47(1, 0) / s6);
        sols[nsol][0] = th5;
        sols[nsol][1] = th6;
        sols[nsol][2] = th7;
        nsol++;
    }
    if (nsol == 0) {
        sols[0][0] = 0.0; sols[0][1] = 0.0; sols[0][2] = 0.0;
        nsol = 1;
    }
    return nsol;
}

// 针对一个给定臂角 theta0 的完整 IK (py: _ik_one_arm_angle)
// 返回解的数量, q_out 按 [q1,q2,q3,q4,q5,q6,q7] 弧度存储
static int ik_one_arm_angle_py(
    const Eigen::Matrix3d& R_des, const Eigen::Vector3d& p_des,
    double theta0, double q_out[7], bool quiet)
{
    Eigen::Vector3d S, W, u_sw;
    double q4_abs;
    compute_swe_from_target_py(R_des, p_des, S, W, q4_abs, u_sw);
    double q4 = q4_abs;

    Eigen::Vector3d E = elbow_from_arm_angle_py(S, W, theta0);

    double sols123[2][3];
    int n123 = solve_q123_from_swe_py(E, W, q4, sols123);
    if (n123 == 0) return 0;

    double q1 = sols123[0][0];
    double q2 = sols123[0][1];
    double q3 = sols123[0][2];

    double q_deg[7] = {0,0,0,0,0,0,0};
    q_deg[0] = q1 * 180.0 / M_PI;
    q_deg[1] = q2 * 180.0 / M_PI;
    q_deg[2] = q3 * 180.0 / M_PI;
    q_deg[3] = q4 * 180.0 / M_PI;

    double T04[4][4];
    forward_kinematics(q_deg, T04, 4);
    Eigen::Matrix3d R_0_4;
    R_0_4 << T04[0][0], T04[0][1], T04[0][2],
             T04[1][0], T04[1][1], T04[1][2],
             T04[2][0], T04[2][1], T04[2][2];

    Eigen::Matrix3d R_4_7 = R_0_4.transpose() * R_des;

    double sols567[2][3];
    int n567 = extract_567_from_T47_py(R_4_7, sols567);

    if (n567 == 0) {
        if (!quiet) printf("ik_one_arm_angle_py: no wrist solution\n");
        return 0;
    }

    double q5 = sols567[0][0];
    double q6 = sols567[0][1];
    double q7 = sols567[0][2];

    q_out[0] = q1;
    q_out[1] = q2;
    q_out[2] = q3;
    q_out[3] = q4;
    q_out[4] = q5;
    q_out[5] = q6;
    q_out[6] = q7;

    return 1;
}

// 完整 IK 入口 (py 风格): 返回最佳解 (弧度)
// 内部暴力扫描 [-pi, pi], step=0.01 找可行臂角
static bool analytical_ik_py(
    const double T_target[4][4],
    const double q_init_deg[7],
    double q_out_deg[7],
    bool quiet = false)
{
    Eigen::Matrix3d R_des;
    Eigen::Vector3d p_des;
    getPoseFromArray(T_target, R_des, p_des);

    double best_cost = 1e30;
    double best_q[7] = {0};
    bool found = false;

    const double step = 0.01;
    int trials = 0;
    for (double theta0 = -M_PI; theta0 < M_PI; theta0 += step) {
        double q_rad[7];
        int ret = ik_one_arm_angle_py(R_des, p_des, theta0, q_rad, true);
        if (ret == 0) continue;
        trials++;

        bool valid = true;
        for (int i = 0; i < 7; i++) {
            if (q_rad[i] < Q_MIN[i] || q_rad[i] > Q_MAX[i]) {
                valid = false;
                break;
            }
        }
        if (!valid) continue;

        double cost = 0.0;
        for (int i = 0; i < 7; i++) {
            double q_init_rad = q_init_deg[i] * M_PI / 180.0;
            double dq = q_rad[i] - q_init_rad;
            cost += dq * dq;
        }
        if (cost < best_cost) {
            best_cost = cost;
            for (int i = 0; i < 7; i++) best_q[i] = q_rad[i];
            found = true;
        }
    }

    if (!found) return false;

    for (int i = 0; i < 7; i++) {
        q_out_deg[i] = best_q[i] * 180.0 / M_PI;
    }
    if (!quiet) {
        printf("py_ik: scanned %d feasible arm angles, best cost=%.6f\n", trials, best_cost);
    }
    return true;
}

// ============================================================
// Py 移植结束
// ============================================================

/**
 * 灏嗚噦瑙掗檺鍒跺埌鏈€杩戠殑鍙鍖洪棿
 * @param 褰撳墠鑷傝锛坮ad锛?
 * @return best_psi 鍚堢悊鐨勮噦瑙掑尯闂达紙rad锛?
 * 
 */
static inline double clampArmAngle(double psi)
{
    // 褰掍竴鍖栧埌 [-蟺, 蟺)
    while (psi < -M_PI) psi += 2.0 * M_PI;
    while (psi >=  M_PI) psi -= 2.0 * M_PI;

    // 宸插湪鍖洪棿1
    if (psi >= ARM_ANGLE_DISTRICT_1_1 && psi <= ARM_ANGLE_DISTRICT_1_2)
        return psi;

    // 宸插湪鍖洪棿2
    if (psi >= ARM_ANGLE_DISTRICT_2_1 && psi <= ARM_ANGLE_DISTRICT_2_2)
        return psi;

    // 涓嶅湪浠讳綍鍖洪棿 鈫?鎵炬渶杩戠殑杈圭晫
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

// 鍙鍏宠妭2锛?锛?妫€鏌ワ紝鍏朵綑鍏宠妭鏉冮噸涓?
static double singularity_penalty(const double q[7]) {
    // 閽熷舰鎯╃綒锛氬湪濂囧紓瑙掑害澶勫嘲鍊间负1锛宻igma鎺у埗瀹藉害
    // penalty = exp(-x^2 / (2*sigma^2))锛宻igma=0.15rad鈮?.6搴?
    const double sigma2 = 0.15 * 0.15;
    const int sing_idx[] = {1, 3, 5};
    const double sing_angles[] = {0.0, M_PI, -M_PI};  // 姣忎釜鍏宠妭鐨勫寮傝

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
 * @param q 璁＄畻鍚庣殑瑙掑害
 * @param q_init 褰撳墠鍒濆瑙掑害
 * @return score 鎵撳垎绯荤粺
 */
static double score_solution(const double q[7],const double q_init[7],
                              const double q_prev[7]) {
    const double w_dq   = 1.0;
    const double w_vel  = 0.5;
    const double w_sing = 5.0;

    double score = 0.0;

    for (int i = 0; i < 7; i++) {
        // 瓒呭嚭闄愪綅鐩存帴娣樻卑
        if (q[i] < Q_MIN[i]*180/M_PI || q[i] > Q_MAX[i]*180/M_PI) return INFINITY;

        // 1. 鍏宠妭瑙掑害鍙樺寲
        double dq = q[i] - q_init[i];
        score += w_dq * dq * dq;

        // 2. 閫熷害鍙樺寲锛堝姞閫熷害浠ｇ悊锛夛細(v_new - v_old) 鍏朵腑 v = dq/dt
        double v_new = (q[i]    - q_init[i]) / DT;  // 褰撳墠姝ラ€熷害
        double v_old = (q_init[i] - q_prev[i]) / DT; // 涓婁竴姝ラ€熷害
        double dv = v_new - v_old;
        score += w_vel * dv * dv;}

    // 3. 濂囧紓鎬э紙閽熷舰锛屽彧鍦ㄥ寮傝搴﹂檮杩戞墠澶э級
    score += w_sing * singularity_penalty(q);

    return score;
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

// 榛勯噾鍒嗗壊 1D 鎼滅储锛屽湪 [lo, hi] 鍐呮壘鏈€灏忓€?
// 杩斿洖鏈€浣庡緱鍒嗭紝psi_best 涓烘渶浼樿噦瑙掞紝q_best 涓哄搴斿叧鑺傝
static double golden_section_1d(double lo, double hi, double tol,
                                const double T_target[4][4],
                                const double q_init[N_JOINTS],
                                const double q_prev[N_JOINTS],
                                double q_best[N_JOINTS],
                                double& psi_best) {
    const double phi   = 0.6180339887498949;   // 1/蠁
    const double phi_c = 1.0 - phi;             // ~0.382

    double a = lo, b = hi;
    double x1 = a + phi_c * (b - a);   // 宸﹀唴鐐?
    double x2 = a + phi   * (b - a);   // 鍙冲唴鐐?

    double q1[N_JOINTS], q2[N_JOINTS];
    double f1 = eval_psi(x1, T_target, q_init, q_prev, q1);
    double f2 = eval_psi(x2, T_target, q_init, q_prev, q2);
    int evals = 2;

    while ((b - a) > tol && evals < 20) {
        if (f1 < f2) {
            // 鎶涘純鍙虫 [x2, b]锛屾棫 x1 鍙樻垚鏂?x2锛堝鐢級
            b  = x2;
            x2 = x1;   f2 = f1;
            x1 = a + phi_c * (b - a);
            f1 = eval_psi(x1, T_target, q_init, q_prev, q1);
        } else {
            // 鎶涘純宸︽ [a, x1]锛屾棫 x2 鍙樻垚鏂?x1锛堝鐢級
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




// q_prev: 涓婁笂鏃跺埢鍏宠妭瑙掞紙鐢ㄤ簬浼拌鍔犻€熷害锛?
/**
 * @param q_init 涓婁竴鏃跺埢锛坱-1锛夌殑鍏宠妭瑙?
 * @param q_prev 涓婁笂鏃跺埢锛坱-2锛夌殑鍏宠妭瑙?
 * @param psi_center 涓婁竴鏃跺埢锛坱-1锛夌殑鑷傝
 * @param q_best 褰撳墠鏃跺埢锛坱锛夌殑鍊欓€夎В
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
        double psi_sample = psi_center + (-7 + i) * step;  // 7鐐瑰绉伴噰鏍?
        double psi = clampArmAngle(psi_sample);//闄愬埗鑷傝
        printf("*******************采样臂角为：%f \n",psi*180/M_PI);
        double q_cand[N_JOINTS];
        analytical_ik_paper(T_target, q_init, psi, q_cand);

        // NaN/Inf 妫€鏌?
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
    double T_target[4][4] = 
        {{-0.02385176519, 0.7667205121, -0.6415378006, 251.5004571},
        {0.05935864847, 0.6416743366, 0.7646767922, 885.2382746},
        {0.9979517244, -0.01984192549, -0.06081655896, 316.4981235},
        {0, 0, 0, 1}};
    double joint_init[7] = {52, 65, 42, 51, 62, 50.5, 85};
    double joint_outcome[7] = {50, 60, 40, 50, 60, 50, 80};

    printf("================================================\n");
    printf("  CPP 原始版 (golden section) VS PY 移植版 (暴力扫描)\n");
    printf("================================================\n\n");

    // ========== C++ 原始版 ==========
    printf("---------- C++ 原始版 ----------\n");
    double joint_cpp[7] = {};
    auto t1 = std::chrono::high_resolution_clock::now();
    int found_cpp = select_optimal_ik_golden(T_target, joint_init, joint_init, joint_cpp);
    auto t2 = std::chrono::high_resolution_clock::now();
    double elapsed_cpp = std::chrono::duration<double, std::micro>(t2 - t1).count();

    for (int i = 0; i < 7; i++) {
        printf("  joint%d = %10.4f deg\n", i, joint_cpp[i]);
    }
    double T_cpp[4][4] = {};
    forward_kinematics(joint_cpp, T_cpp, 7);
    printf("  FK error: ");
    compareTransforms(T_target, T_cpp);
    printf("  time: %.3f us (%.3f ms)\n", elapsed_cpp, elapsed_cpp / 1000.0);

    // ========== PY 移植版 ==========
    printf("\n---------- PY 移植版 ----------\n");
    double joint_py[7] = {};
    auto t3 = std::chrono::high_resolution_clock::now();
    bool found_py = analytical_ik_py(T_target, joint_init, joint_py, false);
    auto t4 = std::chrono::high_resolution_clock::now();
    double elapsed_py = std::chrono::duration<double, std::micro>(t4 - t3).count();

    if (found_py) {
        for (int i = 0; i < 7; i++) {
            printf("  joint%d = %10.4f deg\n", i, joint_py[i]);
        }
        double T_py[4][4] = {};
        forward_kinematics(joint_py, T_py, 7);
        printf("  FK error: ");
        compareTransforms(T_target, T_py);
        printf("  time: %.3f us (%.3f ms)\n", elapsed_py, elapsed_py / 1000.0);
    } else {
        printf("  PY IK FAILED (no solution found)\n");
    }

    // ========== 对比总结 ==========
    printf("\n========================================\n");
    printf("  对比总结:\n");
    printf("  C++版耗时: %.3f us\n", elapsed_cpp);
    printf("  PY版耗时:  %.3f us\n", elapsed_py);
    if (found_py) {
        printf("  速度比 (C++/PY): %.2fx\n", elapsed_py / elapsed_cpp);
        printf("  解差异 (mm):\n");
        double diff_max = 0.0;
        for (int i = 0; i < 7; i++) {
            double d = std::abs(joint_cpp[i] - joint_py[i]);
            if (d > diff_max) diff_max = d;
        }
        printf("  最大关节角偏差: %.6f deg\n", diff_max);
    } else {
        printf("  PY版未找到解，无法对比结果\n");
    }
    printf("========================================\n");

    return 0;
}


