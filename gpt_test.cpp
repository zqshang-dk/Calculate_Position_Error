#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//BDSC系下的地心引力常数
const double mu_bds = 3.986004418e14;  //单位m?/s?
//WGS84坐标系下的地心引力常数
const double mu_gps = 3.986005e14;
//BDSC坐标系下的地球自转角速度
const double OMEGA_e_dot = 7.2921150e-5;   //单位rad/s
//WGS84坐标系下的地球自转角速度
const double OMEGA_e_dot_gps = 7.2921151467e-5;
//圆周率
const double PI = 3.1415926535898;
//GPS参考半长轴
const double Aref = 26559710.0;
//GPS参考升交点赤经变化率
const double OMEGA_REF_dot = -2.6e-9 * 2 * PI;

// ------------------- 星历数据结构 -------------------
struct EphBlock
{
    char sys;           // G / C
    int prn;

    int year, month, day, hour, minute;
    double second;

    double a0, a1, a2;

    double IODE, Crs, Deltan, M0;
    double Cuc, e, Cus, sqrtA;
    double Toe, Cic, OMEGA0, Cis;
    double i0, Crc, omega, OMEGAdot;
    double IDOT, week, L2C, L2P;
    double acc, health, TGD, IODC;

    double transmissionTime;
};

// ------------------- 工具：解析 RINEX 科学计数法 -------------------
double parseLineDouble(const string& s, int start, int length)
{
    string sub = s.substr(start, length);
    // RINEX 科学计数法中使用 D 替换 E
    for (char &c : sub) if (c == 'D' || c == 'd') c = 'E';

    return atof(sub.c_str());
}

// ------------------- 解析一行（头行） -------------------
void parseHeadLine(const string& line, EphBlock &eph)
{
    // 示例格式：
    // G 32 2023 11 17 00 00 00.0 -1.234567D-04 -2.345678D-12 0.000000D+00

    eph.sys = line[0];                     // G / C / E ...
    eph.prn = stoi(line.substr(2, 2));

    eph.year   = stoi(line.substr(5, 4));
    eph.month  = stoi(line.substr(10, 2));
    eph.day    = stoi(line.substr(13, 2));
    eph.hour   = stoi(line.substr(16, 2));
    eph.minute = stoi(line.substr(19, 2));
    eph.second = parseLineDouble(line, 22, 5);

    eph.a0 = parseLineDouble(line,  30, 19);
    eph.a1 = parseLineDouble(line,  49, 19);
    eph.a2 = parseLineDouble(line,  68, 19);
}

// ------------------- 从文件中读取导航星历 -------------------
void readNavFile(const string& filename,
                 vector<EphBlock> &gps,
                 vector<EphBlock> &bds)
{
    ifstream fin(filename);
    if (!fin.is_open())
    {
        cerr << "无法打开文件: " << filename << endl;
        return;
    }

    string line;

    // 跳过头部
    while (getline(fin, line))
    {
        if (line.find("END OF HEADER") != string::npos)
            break;
    }

    // 开始读取每个星历块
    while (getline(fin, line))
    {
        if (line.size() < 2) continue;

        char sys = line[0];
        if (sys != 'G' && sys != 'C')  // 只读取 GPS 和 BDS
            continue;

        EphBlock eph;
        parseHeadLine(line, eph);

        string line2, line3, line4, line5, line6, line7, line8;

        getline(fin, line2);
        getline(fin, line3);
        getline(fin, line4);
        getline(fin, line5);
        getline(fin, line6);
        getline(fin, line7);
        getline(fin, line8);

        // 第二行
        eph.IODE   = parseLineDouble(line2,  4, 19);
        eph.Crs    = parseLineDouble(line2, 23, 19);
        eph.Deltan = parseLineDouble(line2, 42, 19);
        eph.M0     = parseLineDouble(line2, 61, 19);

        // 第三行
        eph.Cuc    = parseLineDouble(line3,  4, 19);
        eph.e       = parseLineDouble(line3, 23, 19);
        eph.Cus    = parseLineDouble(line3, 42, 19);
        eph.sqrtA  = parseLineDouble(line3, 61, 19);

        // 第四行
        eph.Toe    = parseLineDouble(line4, 4, 19);
        eph.Cic    = parseLineDouble(line4, 23, 19);
        eph.OMEGA0  = parseLineDouble(line4, 42, 19);
        eph.Cis    = parseLineDouble(line4, 61, 19);

        // 第五行
        eph.i0        = parseLineDouble(line5, 4, 19);
        eph.Crc       = parseLineDouble(line5, 23, 19);
        eph.omega     = parseLineDouble(line5, 42, 19);
        eph.OMEGAdot  = parseLineDouble(line5, 61, 19);

        // 第六行
        eph.IDOT  = parseLineDouble(line6, 4, 19);
        eph.week  = parseLineDouble(line6, 23, 19);
        eph.L2C   = parseLineDouble(line6, 42, 19);
        eph.L2P   = parseLineDouble(line6, 61, 19);

        // 第七行
        eph.acc     = parseLineDouble(line7,  4, 19);
        eph.health  = parseLineDouble(line7, 23, 19);
        eph.TGD     = parseLineDouble(line7, 42, 19);
        eph.IODC    = parseLineDouble(line7, 61, 19);

        // 第八行
        eph.transmissionTime = parseLineDouble(line8, 4, 19);

        // 分类存放
        if (sys == 'G') gps.push_back(eph);
        if (sys == 'C') bds.push_back(eph);
    }

    fin.close();
}

// //计算BDSC下的卫星位置
// void calBDS(EphBlock BDSdata){
//     //计算长半轴
//     double A = (BDSdata.sqrtA) * (BDSdata.sqrtA);

//     //计算卫星平均角速度
//     double n0 = sqrt(mu_bds / (A * A * A));

//     //计算观测历元到参考历元的时间差
//     double tk = BDSdata.transmissionTime - BDSdata.Toe;//(如果tk大于302400，将tk减去604800；如果tk小于-302400，则将tk加上604800。)

//     //改正平均角速度
//     double n = n0 + BDSdata.Deltan;

//     //计算平近点角
//     double Mk = BDSdata.M0 + n * tk;

//     //迭代计算偏近点角
//     Mk = BDSdata.Ek - BDSdata.e * sin(Ek);

//     //计算真近点角
//     double Vk;
//     sin(Vk) = sqrt(1 - BDSdata.e * BDSdata.e) * sin(Ek) / (1 - BDSdata.e * cos(Ek));
//     cos(Vk) = (cos(Ek) - BDSdata.e) / (1 - BDSdata.e * sin(Ek));

//     //计算纬度幅角
//     double Phi_k = Vk + BDSdata.omega;

//     //纬度幅角改正项
//     double deltauk = BDSdata.Cus * sin(2 * Phi_k) + BDSdata.Cuc * cos(2 * Phi_k);
//     double deltark = BDSdata.Crs * sin(2 * Phi_k) + BDSdata.Crc * cos(2 * Phi_k);
//     double deltaik = BDSdata.Cis * sin(2 * Phi_k) + BDSdata.Cic * cos(2 * Phi_k);

//     //计算改正后的纬度幅角
//     double uk = Phi_k + deltauk;

//     //计算改正后的径向
//     double rk = A * (1 - BDSdata.e * cos(Ek)) + deltark;

//     //计算改正后的轨道倾角
//     double ik = BDSdata.i0 + BDSdata.IDOT * tk + deltaik;

//     //计算卫星在轨道平面内的坐标
//     double xk = rk * cos(uk);
//     double yk = rk * sin(uk);

//     //计算历元升交点经度(ECEF系)
//     double OMEGAk = BDSdata.OMEGA0 + (BDSdata.OMEGAdot - OMEGA_e_dot) * tk - OMEGA_e_dot * BDSdata.Toe;

//     //计算MEO/IGSO卫星在BDSC坐标系中的坐标
//     double Xk = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
//     double Yk = xk * sin(OMEGAk) - yk * cos(ik) * cos(OMEGAk);
//     double Zk = yk * sin(ik);

//     //计算历元升交点经度（惯性系）
//     OMEGAk = BDSdata.OMEGA0 + BDSdata.OMEGAdot * tk - OMEGA_e_dot * BDSdata.Toe;
//     //计算GEO卫星在自定义坐标系中的坐标
//     double XGk = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
//     double YGk = xk * sin(OMEGAk) - yk * cos(ik) * cos(OMEGAk);
//     double ZGk = yk * sin(ik);

//     Matrix3d Rx(phi) = [ (1, 0, 0), (0, cos(phi), sin(phi)), (0, -sin(phi), cos(phi)) ];
//     Matrix3d Rz(phi) = [
//         (cos(phi), sin(phi), 0),
//         (-sin(phi), cos(phi), 0),
//         (0, 0, 1),
//     ];

//     Vector3d(Xk, Yk, Zk) = Rz(OMEGA_e_dot * tk) * Rx(-5 * PI / 180) * Vector3d(XGk, YGk, ZGk);
// }

//计算BDSC下的卫星位置
Vector3d calBDS(const EphBlock& BDSdata){
    //计算长半轴
    double A = (BDSdata.sqrtA) * (BDSdata.sqrtA);

    //计算卫星平均角速度
    double n0 = sqrt(mu_bds / (A * A * A));

    //计算观测历元到参考历元的时间差
    double tk = BDSdata.transmissionTime - BDSdata.Toe;//(如果tk大于302400，将tk减去604800；如果tk小于-302400，则将tk加上604800。)
    if(tk>302400.0){
        tk -= 604800.0;
    }
    if(tk<302400.0){
        tk += 604800.0;
    }

    //改正平均角速度
    double n = n0 + BDSdata.Deltan;

    double Mk = BDSdata.M0 + n * tk;

    //计算平近点角迭代求Ek
    double Ek = Mk;
    for (int i = 0; i < 10;i++){
        Ek = Ek + (Mk - Ek + BDSdata.e * sin(Ek));
    }


    //计算真近点角
    double sinVk = sqrt(1 - BDSdata.e * BDSdata.e) * sin(Ek) / (1 - BDSdata.e * cos(Ek));
    double cosVk = (cos(Ek) - BDSdata.e) / (1 - BDSdata.e * sin(Ek));
    double Vk = atan2(sinVk, cosVk);

    //计算纬度幅角
    double Phi_k = Vk + BDSdata.omega;

    //纬度幅角改正项
    double deltauk = BDSdata.Cus * sin(2 * Phi_k) + BDSdata.Cuc * cos(2 * Phi_k);
    double deltark = BDSdata.Crs * sin(2 * Phi_k) + BDSdata.Crc * cos(2 * Phi_k);
    double deltaik = BDSdata.Cis * sin(2 * Phi_k) + BDSdata.Cic * cos(2 * Phi_k);

    //计算改正后的纬度幅角
    double uk = Phi_k + deltauk;

    //计算改正后的径向
    double rk = A * (1 - BDSdata.e * cos(Ek)) + deltark;

    //计算改正后的轨道倾角
    double ik = BDSdata.i0 + BDSdata.IDOT * tk + deltaik;

    //计算卫星在轨道平面内的坐标
    double xk = rk * cos(uk);
    double yk = rk * sin(uk);

    //计算历元升交点经度(ECEF系)
    //double OMEGAk = BDSdata.OMEGA0 + (BDSdata.OMEGAdot - OMEGA_e_dot) * tk - OMEGA_e_dot * BDSdata.Toe;

    double Xk, Yk, Zk;
    if(BDSdata.prn<=5){  //GPS
        double OMEGAk = BDSdata.OMEGA0 + BDSdata.OMEGAdot * tk - OMEGA_e_dot * BDSdata.Toe;
        double XGk = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
        double YGk = xk * sin(OMEGAk) - yk * cos(ik) * cos(OMEGAk);
        double ZGk = yk * sin(ik);

        double phi = -5.0 * PI / 180;
        Matrix3d Rx;
        Rx << 1, 0, 0,
            0, cos(phi), sin(phi),
            0, -sin(phi), cos(phi);

        double psi = OMEGA_e_dot * tk;
        Matrix3d Rz;
        Rz << cos(psi), sin(psi), 0,
            -sin(psi), cos(psi), 0,
            0, 0, 1;

        Vector3d pos = Rz * Rx * Vector3d(XGk, YGk, ZGk);

        Xk = pos(0);
        Yk = pos(1);
        Zk = pos(2);
    }

    else{//MEO/IGSO
        double OMEGAk = BDSdata.OMEGA0 + (BDSdata.OMEGAdot - OMEGA_e_dot) * tk - OMEGA_e_dot * BDSdata.Toe;
        Xk = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
        Yk = xk * sin(OMEGAk) - yk * cos(ik) * cos(OMEGAk);
        Zk = yk * sin(ik);

    }
    return Vector3d(Xk, Yk, Zk);
}


//计算GPS卫星的位置
Vector3d calGPS(const EphBlock& GPSdata){
    // ---------------------- 步骤1：计算时间差tk ----------------------
    double tk = GPSdata.transmissionTime - GPSdata.Toe;
    // 处理周跳（若tk超出±302400秒，调整±604800秒）
    if (tk > 302400.0)
        tk -= 604800.0;
    if (tk < -302400.0)
        tk += 604800.0;

    // ---------------------- 步骤2：计算轨道半长轴及平均角速度 ----------------------
    double A = GPSdata.sqrtA * GPSdata.sqrtA; // 半长轴 A = (sqrtA)^2
    double n0 = sqrt(mu_gps / (A * A * A));       // 计算平均角速度n0
    double n = n0 + GPSdata.Deltan;           // 修正后的平均角速度


    // ---------------------- 步骤3：解开普勒方程（计算偏近点角Ek） ----------------------
    double Mk = GPSdata.M0 + n * tk; // 平近点角
    double Ek = Mk;                  // 初始值E0 = Mk
    // 迭代求解（至少3次）
    for (int i = 0; i < 10; ++i) {
        Ek = Ek + (Mk - Ek + GPSdata.e * sin(Ek)) / (1 - GPSdata.e * cos(Ek));
    }


    // ---------------------- 步骤4：计算真近点角vk ----------------------
    double sqrt_1me2 = sqrt(1 - GPSdata.e * GPSdata.e);
    double vk = 2 * atan2( sqrt_1me2 * sin(Ek), (cos(Ek) - GPSdata.e) );


    // ---------------------- 步骤5：计算纬度幅角及二阶调和摄动 ----------------------
    double Phik = vk + GPSdata.omega; // 纬度幅角

    // 二阶调和摄动修正项
    double sin2Phik = sin(2 * Phik);
    double cos2Phik = cos(2 * Phik);
    double deltaUk = GPSdata.Cus * sin2Phik + GPSdata.Cuc * cos2Phik; // 纬度幅角修正
    double deltaRk = GPSdata.Crs * sin2Phik + GPSdata.Crc * cos2Phik; // 径向修正
    double deltaIk = GPSdata.Cis * sin2Phik + GPSdata.Cic * cos2Phik; // 轨道倾角修正


    // ---------------------- 步骤6：修正轨道参数 ----------------------
    double uk = Phik + deltaUk;      // 修正后纬度幅角
    double rk = A * (1 - GPSdata.e * cos(Ek)) + deltaRk; // 修正后轨道半径
    double ik = GPSdata.i0 + GPSdata.IDOT * tk + deltaIk; // 修正后轨道倾角


    // ---------------------- 步骤7：计算升交点赤经Ωk ----------------------
    double OmegaDot = OMEGA_REF_dot + GPSdata.OMEGAdot; // 升交点赤经变化率
    double OMEGAk = GPSdata.OMEGA0 + (OmegaDot - OMEGA_e_dot_gps) * tk - OMEGA_e_dot_gps * GPSdata.Toe;


    // ---------------------- 步骤8：轨道平面内坐标(xk', yk') ----------------------
    double xk_prime = rk * cos(uk);
    double yk_prime = rk * sin(uk);


    // ---------------------- 步骤9：地球固连坐标系（ECEF）下的卫星位置 ----------------------
    double cosOmega = cos(OMEGAk);
    double sinOmega = sin(OMEGAk);
    double cosi = cos(ik);
    double sini = sin(ik);

    double xk = xk_prime * cosOmega - yk_prime * cosi * sinOmega;
    double yk = xk_prime * sinOmega + yk_prime * cosi * cosOmega;
    double zk = yk_prime * sini;

    return Vector3d(xk, yk, zk);
}





// // ------------------- 输出到 TXT -------------------
// void saveTxt(const string& filename, const vector<EphBlock>& data)
// {
//     ofstream fout(filename);

//     for (const auto &e : data)
//     {
//         fout << fixed << setprecision(10);
//         fout << e.sys << setw(2) << e.prn << " "
//              << e.year << "-" << e.month << "-" << e.day << " "
//              << e.hour << ":" << e.minute << ":" << e.second << "\n";

//         fout << "a0 = " << e.a0 << "  a1 = " << e.a1 << "  a2 = " << e.a2 << "\n";
//         fout << "e = " << e.e << "  sqrtA = " << e.sqrtA << "  Toe = " << e.Toe << "\n";
//         fout << "i0 = " << e.i0 << "  omega = " << e.omega << "\n";
//         fout << "---------------------------------------\n";
//     }

//     fout.close();
// }
// ------------------- 保存 GPS 位置到文件 -------------------
void saveGPSPosition(const string& filename, const vector<EphBlock>& gps)
{
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "无法创建 GPS 位置文件: " << filename << endl;
        return;
    }
    fout << fixed << setprecision(3);  // 位置一般保留到毫米级（3位小数，单位：米）
    fout << "# PRN        X (m)              Y (m)              Z (m)" << endl;

    for (const auto& g : gps) {
        Vector3d pos = calGPS(g);
        fout << "G" << setw(2) << g.prn << "   "
             << setw(14) << pos(0) << "   "
             << setw(14) << pos(1) << "   "
             << setw(14) << pos(2) << endl;
    }
    fout.close();
    cout << "GPS 卫星位置已保存至: " << filename << endl;
}

// ------------------- 保存 BDS 位置到文件 -------------------
void saveBDSPosition(const string& filename, const vector<EphBlock>& bds)
{
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "无法创建 BDS 位置文件: " << filename << endl;
        return;
    }
    fout << fixed << setprecision(3);
    fout << "# PRN        X (m)              Y (m)              Z (m)" << endl;

    for (const auto& b : bds) {
        Vector3d pos = calBDS(b);
        fout << "C" << setw(2) << b.prn << "   "
             << setw(14) << pos(0) << "   "
             << setw(14) << pos(1) << "   "
             << setw(14) << pos(2) << endl;
    }
    fout.close();
    cout << "BDS 卫星位置已保存至: " << filename << endl;
}

// ------------------- 主函数 -------------------
int main()
{
    string filename = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/brdm3350.19p";  // 你可以替换成任意 RINEX 文件

    vector<EphBlock> gps, bds;

    readNavFile(filename, gps, bds);

    cout << "开始计算卫星位置并保存..." << endl;

    // 保存 GPS 位置（同时会在控制台看到进度）
    saveGPSPosition("gps_position.txt", gps);

    // 保存 BDS 位置
    saveBDSPosition("bds_position.txt", bds);

    // saveTxt("gps_ephemeris.txt", gps);
    // saveTxt("bds_ephemeris.txt", bds);

    cout << "读取完成！" << endl;
    cout << "GPS 星历数目: " << gps.size() << endl;
    cout << "BDS 星历数目: " << bds.size() << endl;

    return 0;
}
