#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include <set>
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
const double C_LIGHT = 299792458;

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
    eph.prn = stoi(line.substr(1, 3));

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

// 将公历日期转换为GPS周和周内秒
// 输入: year, month, day, hour, minute, second
// 输出: gps_week (int, GPS周号), week_second (double, 周内秒)
// 返回0成功, -1错误
int CalendarToGpsTime(int year, int month, int day, int hour, int minute, double second, int& gps_week, double& week_second) {
    if (year < 1980 || month < 1 || month > 12 || day < 1 || day > 31 || hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0 || second >= 60) {
        return -1;  // 无效日期
    }

    // GPS起始历元: 1980-01-06 00:00:00 UTC
    // 计算从1980-01-06到目标日期的总天数
    int total_days = 0;

    // 计算完整年份天数 (1980 to year-1)
    for (int y = 1980; y < year; ++y) {
        total_days += (y % 4 == 0 && (y % 100 != 0 || y % 400 == 0)) ? 366 : 365;
    }

    // 当前年月份天数
    int month_days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)) month_days[1] = 29;  // 闰年

    for (int m = 1; m < month; ++m) {
        total_days += month_days[m - 1];
    }
    total_days += day - 6;  // 从1980-01-06开始, 减去前5天 (Jan 1-5)

    // 如果total_days <0, 错误
    if (total_days < 0) return -1;

    // GPS周号和周内天
    gps_week = total_days / 7;
    int day_of_week = total_days % 7;

    // 周内秒
    week_second = day_of_week * 86400.0 + hour * 3600.0 + minute * 60.0 + second;

    return 0;  // 成功
}

// 现在的钟差计算还是错误的，钟差计算有两个需要注意的点
// （1）最后一定要转化为时间
// （2）toc是指卫星参考时刻，
// 从年月日转换为GPS周和周内秒

Vector4d calBDS(const EphBlock& BDSdata,double t){
    //计算长半轴
    double A = (BDSdata.sqrtA) * (BDSdata.sqrtA);

    //计算卫星平均角速度
    double n0 = sqrt(mu_bds / (A * A * A));

    //计算观测历元到参考历元的时间差
    double tk = t - BDSdata.Toe;//(如果tk大于302400，将tk减去604800；如果tk小于-302400，则将tk加上604800。)
    if(tk>302400.0){
        tk -= 604800.0;
    }
    if(tk<-302400.0){
        tk += 604800.0;
    }

    //计算钟差
    double dT = (BDSdata.a0 + BDSdata.a1 * tk + BDSdata.a2 * tk * tk);
    //改正平均角速度
    double n = n0 + BDSdata.Deltan;

    double Mk = BDSdata.M0 + n * tk;

    //计算平近点角迭代求Ek
    double Ek = Mk;
    for (int i = 0; i < 10;i++){
        //Ek = Ek + (Mk - Ek + BDSdata.e * sin(Ek));
        Ek = Mk + BDSdata.e * sin(Ek);
    }


    //计算真近点角
    double sinVk = sqrt(1 - BDSdata.e * BDSdata.e) * sin(Ek) / (1 - BDSdata.e * cos(Ek));
    double cosVk = (cos(Ek) - BDSdata.e) / (1 - BDSdata.e * cos(Ek));
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

    double X_Inertial, Y_Inertial, Z_Inertial; // 惯性系下的坐标
    double Xk, Yk, Zk; // 最终 ECEF 坐标

    double psi = OMEGA_e_dot * tk;
    Matrix3d Rz;
    Rz << cos(psi), sin(psi), 0,
        -sin(psi), cos(psi), 0,
        0, 0, 1;
    if(BDSdata.prn<=5){  //GEO卫星
        double OMEGAk = BDSdata.OMEGA0 + BDSdata.OMEGAdot * tk - OMEGA_e_dot * BDSdata.Toe;
        //double OMEGAk = BDSdata.OMEGA0 + (BDSdata.OMEGAdot - OMEGA_e_dot) * tk - OMEGA_e_dot * BDSdata.Toe;
        X_Inertial = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
        Y_Inertial = xk * sin(OMEGAk) + yk * cos(ik) * cos(OMEGAk);
        Z_Inertial = yk * sin(ik);

        double phi = -5.0 * PI / 180;
        Matrix3d Rx;
        Rx << 1, 0, 0,
            0, cos(phi), sin(phi),
            0, -sin(phi), cos(phi);

        Vector3d pos = Rz * Rx * Vector3d(X_Inertial, Y_Inertial, Z_Inertial);

        Xk = pos(0);
        Yk = pos(1);
        Zk = pos(2);
    }

    else{//MEO/IGSO
        double OMEGAk = BDSdata.OMEGA0 + (BDSdata.OMEGAdot - OMEGA_e_dot) * tk - OMEGA_e_dot * BDSdata.Toe;
        X_Inertial = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
        Y_Inertial = xk * sin(OMEGAk) + yk * cos(ik) * cos(OMEGAk);
        Z_Inertial = yk * sin(ik);

        Vector3d pos_ECEF = Rz.transpose() * Vector3d(X_Inertial, Y_Inertial, Z_Inertial);
        Xk = pos_ECEF(0);
        Yk = pos_ECEF(1);
        Zk = pos_ECEF(2);
    }
    return Vector4d(Xk, Yk, Zk,dT);
}


Vector4d calGPS(const EphBlock& GPSdata,double t){
    double tk = t - GPSdata.Toe;
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;

    //计算钟差
    //double deltat=GPSdata.transmissionTime-
    double dT = (GPSdata.a0 + GPSdata.a1 * tk + GPSdata.a2 * tk * tk);

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
    double vk = atan2( sqrt_1me2 * sin(Ek), (cos(Ek) - GPSdata.e) );


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
    double OMEGAk = GPSdata.OMEGA0 + (GPSdata.OMEGAdot - OMEGA_e_dot_gps) * tk - OMEGA_e_dot_gps * GPSdata.Toe;

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

    // 添加相对论修正 (只对GPS, 单位秒)
    double F = -4.442807633e-10;  // -2*sqrt(mu_gps)/c
    double rel_corr = F * GPSdata.e * GPSdata.sqrtA * sin(Ek);
    dT += rel_corr;

    return Vector4d(xk, yk, zk,dT);
}

// 为给定PRN和t, 找最佳EphBlock (closest Toe, |tk| < 7200s)
const EphBlock* findBestEph(const vector<EphBlock>& sats, int prn, double t) {
    const EphBlock* best = nullptr;
    double min_diff = 7200.0;  // 最大有效范围 (2小时)

    for (const auto& s : sats) {
        if (s.prn == prn) {
            double diff = fabs(t - s.Toe);
            if (diff < min_diff) {
                min_diff = diff;
                best = &s;
            }
        }
    }
    return best;  // 如果无, 返回nullptr
}
// ------------------- 主函数 -------------------
int main()
{
    string filename = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/brdm3350.19p";  

    vector<EphBlock> gps, bds;

    readNavFile(filename, gps, bds);
    // 假设文件日期2019-12-01, 计算GPS周
    int gps_week;
    double dummy_second;  // 不用
    if (CalendarToGpsTime(2019, 12, 1, 0, 0, 0.0, gps_week, dummy_second) != 0) {
        cerr << "日期转换错误!" << endl;
        return 1;
    }

    // 收集独特PRN
    set<int> gps_prns, bds_prns;
    for (const auto& g : gps) gps_prns.insert(g.prn);
    for (const auto& b : bds) bds_prns.insert(b.prn);

    // 输出文件
    ofstream gps_out("gps_position.txt");
    ofstream bds_out("bds_position.txt");

    gps_out << fixed << setprecision(3);
    bds_out << fixed << setprecision(3);

    // 每15min (900s), 从t=0到86400
    double interval = 900.0;
    for (double t = 0.0; t < 86400.0; t += interval) {
        // 计算当前UTC时间 (year month day hour min sec)
        int hour = static_cast<int>(t / 3600);
        int min = static_cast<int>((t - hour*3600) / 60);
        double sec = t - hour*3600 - min*60;

        // 输出* 行
        gps_out << "* 2019 12 1 " << hour << " " << min << " " << fixed << setprecision(8) << sec << " " << interval << endl;
        bds_out << "* 2019 12 1 " << hour << " " << min << " " << fixed << setprecision(8) << sec << " " << interval << endl;

        // GPS
        for (int prn : gps_prns) {
            const EphBlock* eph = findBestEph(gps, prn, t);
            if (eph) {
                Vector4d pos = calGPS(*eph, t);
                gps_out << "G" << setfill('0') << setw(2) << prn << " "
                        << setw(13) << pos(0) << " " << setw(13) << pos(1) << " " << setw(13) << pos(2) << " "
                        << setprecision(8) << pos(3) << endl;
            }
        }

        // BDS
        for (int prn : bds_prns) {
            const EphBlock* eph = findBestEph(bds, prn, t);
            if (eph) {
                Vector4d pos = calBDS(*eph, t);
                bds_out << "C" << setfill('0') << setw(2) << prn << " "
                        << setw(13) << pos(0) << " " << setw(13) << pos(1) << " " << setw(13) << pos(2) << " "
                        << setprecision(8) << pos(3) << endl;
            }
        }
    }

    gps_out.close();
    bds_out.close();

        // ==================== 新增部分：生成每60秒的轨迹文件（用于绘图） ====================

    cout << "正在生成全天轨迹文件（每60秒一点）用于绘图..." << endl;

    // 打开CSV文件
    ofstream gps_track("gps_track.csv");
    ofstream bds_track("bds_track.csv");

    // 写表头
    gps_track << "PRN,t,X,Y,Z\n";
    bds_track << "PRN,t,X,Y,Z\n";

    gps_track << fixed << setprecision(6);
    bds_track << fixed << setprecision(6);

    // 每60秒计算一次（全天 1440 个点）
    for (double t = 0.0; t < 86400.0; t += 60.0) {
        // GPS 所有卫星
        for (int prn : gps_prns) {
            const EphBlock* eph = findBestEph(gps, prn, t);
            if (eph) {
                Vector4d pos = calGPS(*eph, t);
                gps_track << prn << ","
                          << t << ","
                          << pos(0) << ","
                          << pos(1) << ","
                          << pos(2) << "\n";
            }
        }

        // BDS 所有卫星
        for (int prn : bds_prns) {
            const EphBlock* eph = findBestEph(bds, prn, t);
            if (eph) {
                Vector4d pos = calBDS(*eph, t);
                bds_track << prn << ","
                          << t << ","
                          << pos(0) << ","
                          << pos(1) << ","
                          << pos(2) << "\n";
            }
        }

        // 进度提示（每100个点打印一次）
        if (static_cast<int>(t) % 6000 == 0) {  // 每10分钟提示一次
            int hour = static_cast<int>(t / 3600);
            int minute = static_cast<int>((t - hour*3600)/60);
            cout << "   已计算到: " << hour << ":" << setw(2) << setfill('0') << minute << endl;
        }
    }

    gps_track.close();
    bds_track.close();

    cout << "轨迹文件生成完成！" << endl;
    cout << "   gps_track.csv  → 用于 GPS 轨迹绘图（每60秒）" << endl;
    cout << "   bds_track.csv  → 用于 BDS 轨迹绘图（每60秒）" << endl;
    // =====================================================================

    cout << "输出完成! GPS文件: gps_position.txt, BDS文件: bds_position.txt" << endl;
    cout << "GPS PRN数: " << gps_prns.size() << ", BDS PRN数: " << bds_prns.size() << endl;

    return 0;

}
