#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <iomanip>
#include <windows.h>

using namespace std;

// 常数定义
const double PI = 3.14159265358979323846;
const double MU_GPS = 3.986005e14;        // WGS84地球引力常数 (m^3/s^2)
const double MU_BDS = 3.986004418e14;    // BDCS地球引力常数 (m^3/s^2)
const double OMEGA_E_GPS = 7.2921151467e-5;  // WGS84地球自转角速度 (rad/s)
const double OMEGA_E_BDS = 7.2921150e-5;     // BDCS地球自转角速度 (rad/s)

// 基本广播星历块结构体
struct EPHEMERISBLOCK {
    // PRN号和卫星系统
    int PRN;
    char satSystem;  // 'G' for GPS, 'C' for BDS
    
    // 星历参考时间
    int year, month, day, hour, minute;
    double second;
    double toc;  // 钟差参考时间 (GPS秒)
    double toe;  // 星历参考时间 (GPS秒)
    
    // 钟差参数
    double a0, a1, a2;
    
    // 轨道参数
    double IODE, Crs, Deltan, M0;        // ORBIT - 1
    double Cuc, e, Cus, SqrtA;           // ORBIT - 2  
    double Cic, OMEGA, Cis;              // ORBIT - 3
    double i0, Crc, omega, OMEGAdot;     // ORBIT - 4
    double IDOT, GpsWeekNumber, L2C, L2P; // ORBIT - 5
    double SatAccuracy, SatHealth, TGD, IODC; // ORBIT - 6
};

struct SatellitePosition {
    int PRN;
    char satSystem;
    double X, Y, Z;  // ECEF坐标 (m)
    double time;      // 时间 (GPS秒)
    int hour, minute; // 计算时间的小时和分钟
    double tow;       // GPS周内秒
};

// 精密星历点结构体
struct PreciseEphemerisPoint {
    int PRN;
    char satSystem;
    double X, Y, Z;  // ECEF坐标 (m)
    double time;      // 时间 (GPS秒)
    int hour, minute; // 时间的小时和分钟
};

// 位置比较结果结构体
struct PositionComparison {
    int PRN;
    char satSystem;
    double time;           // 时间 (GPS秒)
    double broadcastX, broadcastY, broadcastZ;
    double preciseX, preciseY, preciseZ;
    double dX, dY, dZ;    // 差值 (广播-精密)
    double distance;      // 位置差模长
    int hour, minute;     // 时间
};

// 从年月日转换为GPS周和周内秒
int Calendar2GpsTime(int nYear, int nMonth, int nDay, int nHour, int nMinute, double dSecond, double &WeekSecond) {
    if (nYear < 1980 || nMonth < 1 || nMonth > 12 || nDay < 1 || nDay > 31) {
        return -1;
    }
    
    // 计算从1980年1月6日到当前日期的天数
    int totalDays = 0;
    
    // 计算从1980年到当前年前一年的总天数
    for (int year = 1980; year < nYear; year++) {
        if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
            totalDays += 366;
        } else {
            totalDays += 365;
        }
    }
    
    // 计算当前年的天数
    int monthDays[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if ((nYear % 4 == 0 && nYear % 100 != 0) || (nYear % 400 == 0)) {
        monthDays[1] = 29;
    }
    
    for (int month = 1; month < nMonth; month++) {
        totalDays += monthDays[month - 1];
    }
    
    totalDays += (nDay - 6);  // GPS时间从1980年1月6日开始
    
    int weekno = totalDays / 7;
    int dayofweek = totalDays % 7;
    
    WeekSecond = dayofweek * 86400.0 + nHour * 3600.0 + nMinute * 60.0 + dSecond;
    
    return weekno;
}

// 读取广播星历文件（只读取GPS和BDS）
int ReadBrodcastEphemeris(const char* filename, vector<EPHEMERISBLOCK>& ephemerisBlocks) {
    FILE* pfEph = fopen(filename, "r");
    if (!pfEph) {
        cerr << "Error: Cannot open file " << filename << endl;
        return -1;
    }
    
    char strLine[256];
    
    // 跳过文件头
    while (fgets(strLine, sizeof(strLine), pfEph)) {
        if (strstr(strLine, "END OF HEADER")) {
            break;
        }
    }
    
    int totalRead = 0;
    int gpsRead = 0;
    int bdsRead = 0;
    
    while (fgets(strLine, sizeof(strLine), pfEph)) {
        char satType;
        int prn, year, month, day, hour, minute;
        double second, a0, a1, a2;
        
        if (sscanf(strLine, "%c%2d %d %d %d %d %d %lf %lf %lf %lf", 
                   &satType, &prn, &year, &month, &day, &hour, &minute, &second, &a0, &a1, &a2) >= 9) {
            
            // 只处理GPS和BDS卫星
            if (satType != 'G' && satType != 'C') {
                // 跳过后续7行
                for (int i = 0; i < 7; i++) {
                    fgets(strLine, sizeof(strLine), pfEph);
                }
                continue;
            }
            
            totalRead++;
            
            EPHEMERISBLOCK eph;
            eph.satSystem = satType;
            eph.PRN = prn;
            eph.year = year;
            eph.month = month;
            eph.day = day;
            eph.hour = hour;
            eph.minute = minute;
            eph.second = second;
            eph.a0 = a0;
            eph.a1 = a1;
            eph.a2 = a2;
            
            // 计算toc (钟差参考时间)
            Calendar2GpsTime(year, month, day, hour, minute, second, eph.toc);
            
            // 读取轨道参数
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &eph.IODE, &eph.Crs, &eph.Deltan, &eph.M0);
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &eph.Cuc, &eph.e, &eph.Cus, &eph.SqrtA);
            
            fgets(strLine, sizeof(strLine), pfEph);
            double toe;
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &toe, &eph.Cic, &eph.OMEGA, &eph.Cis);
            
            // // 计算toe (星历参考时间)
            // eph.toe = eph.toc + (toe - (int)eph.toc % 86400);
            // if (eph.toe - eph.toc > 302400) eph.toe -= 604800;
            // if (eph.toe - eph.toc < -302400) eph.toe += 604800;

            // BDS特殊处理：toe计算
           if (satType == 'C') {
                // BDS卫星：toe需要特殊处理
                // BDS的toe是相对于参考时间的偏移，需要正确计算
                eph.toe = eph.toc + toe;
                
                // 确保toe在合理范围内（±3.5天）
                double timeDiff = eph.toe - eph.toc;
                if (timeDiff > 302400.0) eph.toe -= 604800.0;
                else if (timeDiff < -302400.0) eph.toe += 604800.0;
            } else {
                // GPS卫星计算保持不变
                eph.toe = eph.toc + (toe - fmod(eph.toc, 86400.0));
                if (eph.toe - eph.toc > 302400.0) eph.toe -= 604800.0;
                if (eph.toe - eph.toc < -302400.0) eph.toe += 604800.0;
            }
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &eph.i0, &eph.Crc, &eph.omega, &eph.OMEGAdot);
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &eph.IDOT, &eph.L2C, &eph.GpsWeekNumber, &eph.L2P);
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &eph.SatAccuracy, &eph.SatHealth, &eph.TGD, &eph.IODC);
            
            // 跳过最后一行
            fgets(strLine, sizeof(strLine), pfEph);
            
            ephemerisBlocks.push_back(eph);
            
            if (satType == 'G') gpsRead++;
            else if (satType == 'C') bdsRead++;
        }
    }
    
    fclose(pfEph);
    cout << "Read " << totalRead << " ephemeris blocks (GPS: " << gpsRead << ", BDS: " << bdsRead << ")" << endl;
    return 0;
}

// 读取精密星历文件 (SP3格式，只读取GPS和BDS)
vector<PreciseEphemerisPoint> ReadPreciseEphemeris(const char* filename) {
    vector<PreciseEphemerisPoint> precisePoints;
    
    ifstream sp3file(filename);
    if (!sp3file.is_open()) {
        cerr << "Error: Cannot open precise ephemeris file " << filename << endl;
        return precisePoints;
    }
    
    string line;
    int year, month, day, hour, minute;
    double second;
    double gpsWeekSecond = 0;
    
    // 读取文件头
    while (getline(sp3file, line)) {
        // 检查是否到文件头结束
        if (line.find("%c") == 0 || line.find("%f") == 0 || 
            line.find("%i") == 0 || line.find("/*") == 0) {
            continue;
        }
        
        // 跳过注释行
        if (line[0] == '*') {
            // 解析时间行: *  2019 12  1  0  0  0.00000000
            if (sscanf(line.c_str(), "* %d %d %d %d %d %lf", 
                       &year, &month, &day, &hour, &minute, &second) == 6) {
                // 转换为GPS秒
                Calendar2GpsTime(year, month, day, hour, minute, second, gpsWeekSecond);
            }
            continue;
        }
        
        // 跳过非数据行
        if (line.length() < 4 || line[0] != 'P') {
            continue;
        }
        
        // 解析卫星位置行: PG01  23977.126797  10034.346106  11299.196542
        char satType;
        int prn;
        double x, y, z;
        
        if (sscanf(line.c_str(), "P%c%2d %lf %lf %lf", 
                   &satType, &prn, &x, &y, &z) == 5) {
            // 只处理GPS和BDS卫星
            if (satType != 'G' && satType != 'C') {
                continue;
            }
            
            PreciseEphemerisPoint point;
            point.PRN = prn;
            point.satSystem = satType;
            point.X = x * 1000.0;  // km转换为m
            point.Y = y * 1000.0;
            point.Z = z * 1000.0;
            // point.time = gpsWeekSecond;
            // point.hour = hour;
            // point.minute = minute;
            // BDS时间调整：如果精密星历使用GPS时间，BDS需要调整
            if (satType == 'C') {
                point.time = gpsWeekSecond - 14.0;  // BDT = GPST - 14s
            } else {
                point.time = gpsWeekSecond;
            }

            point.hour = hour;
            point.minute = minute;
            
            precisePoints.push_back(point);
        }
    }
    
    sp3file.close();
    
    // 统计GPS和BDS卫星数量
    int gpsCount = 0, bdsCount = 0;
    for (const auto& point : precisePoints) {
        if (point.satSystem == 'G') gpsCount++;
        else if (point.satSystem == 'C') bdsCount++;
    }
    
    cout << "Read " << precisePoints.size() << " precise ephemeris points (GPS: " 
         << gpsCount << ", BDS: " << bdsCount << ")" << endl;
    
    return precisePoints;
}

// 查找给定时间最近的广播星历
const EPHEMERISBLOCK* FindClosestEphemeris(int prn, char satSystem, double time, 
                                         const vector<EPHEMERISBLOCK>& ephemerisBlocks) {
    const EPHEMERISBLOCK* closestEph = nullptr;
    double minTimeDiff = 1e9;
    
    for (const auto& eph : ephemerisBlocks) {
        if (eph.PRN == prn && eph.satSystem == satSystem) {
            double timeDiff = fabs(eph.toe - time);
            if (timeDiff < minTimeDiff && timeDiff < 7200.0) { // 限制在2小时内
                minTimeDiff = timeDiff;
                closestEph = &eph;
            }
        }
    }
    
    return closestEph;
}

// // GPS/BDS卫星位置计算
// SatellitePosition CalculateSatellitePosition(const EPHEMERISBLOCK& eph, double transmitTime, int hour, int minute) {
//     SatellitePosition pos;
//     pos.PRN = eph.PRN;
//     pos.satSystem = eph.satSystem;
//     pos.time = transmitTime;
//     pos.hour = hour;
//     pos.minute = minute;
//     pos.tow = fmod(transmitTime, 86400.0);
    
//     bool isBDS = (eph.satSystem == 'C');

//     // BDS时间系统调整：BDT = GPST - 14s
//     double bdt_offset = isBDS ? 14.0 : 0.0;
//     double bdtTime = transmitTime - bdt_offset;

//     double mu = isBDS ? MU_BDS : MU_GPS;
//     double omega_e = isBDS ? OMEGA_E_BDS : OMEGA_E_GPS;
    
//     // // BDS时间系统调整：BDS使用BDT，与GPS时间有14秒偏差
//     // double timeOffset = isBDS ? 14.0 : 0.0;
//     // double adjustedTransmitTime = transmitTime - timeOffset;
    
    
//     // 计算时间差
//     double tk = transmitTime - eph.toe;
    
//     // 处理周跳
//     if (tk > 302400.0) tk -= 604800.0;
//     if (tk < -302400.0) tk += 604800.0;

//     // BDS地球自转角速度需要调整（BDCS系）
//     if (isBDS) {
//         omega_e = 7.2921150e-5;  // 确保使用正确的BDS地球自转角速度
//     }

//     // 计算平均角速度
//     double A = eph.SqrtA * eph.SqrtA;
//     double n0 = sqrt(mu / (A * A * A));
//     double n = n0 + eph.Deltan;
    
//     // 计算平近点角
//     double Mk = eph.M0 + n * tk;
    
//     // 迭代计算偏近点角
//     double Ek = Mk;
//     double Ek_prev;
//     for (int i = 0; i < 10; i++) {
//         Ek_prev = Ek;
//         Ek = Mk + eph.e * sin(Ek);
//         if (fabs(Ek - Ek_prev) < 1e-12) break;
//     }
    
//     // 计算真近点角
//     double sin_vk = (sqrt(1 - eph.e * eph.e) * sin(Ek)) / (1 - eph.e * cos(Ek));
//     double cos_vk = (cos(Ek) - eph.e) / (1 - eph.e * cos(Ek));
//     double vk = atan2(sin_vk, cos_vk);
    
//     // 计算纬度幅角
//     double phik = vk + eph.omega;
    
//     // 计算周期改正项
//     double delta_uk = eph.Cus * sin(2 * phik) + eph.Cuc * cos(2 * phik);
//     double delta_rk = eph.Crs * sin(2 * phik) + eph.Crc * cos(2 * phik);
//     double delta_ik = eph.Cis * sin(2 * phik) + eph.Cic * cos(2 * phik);
    
//     // 计算改正后的参数
//     double uk = phik + delta_uk;
//     double rk = A * (1 - eph.e * cos(Ek)) + delta_rk;
//     double ik = eph.i0 + eph.IDOT * tk + delta_ik;
    
//     // 计算轨道平面坐标
//     double xk = rk * cos(uk);
//     double yk = rk * sin(uk);
    
//     // // 计算升交点经度
//     // double OMEGAk = eph.OMEGA + (eph.OMEGAdot - omega_e) * tk - omega_e * eph.toe;
    
//     // BDS坐标系处理：GEO卫星有特殊处理，但这里先处理MEO/IGSO
//     double OMEGAk;
//     if (isBDS) {
//         // BDS使用BDCS坐标系，需要特殊处理
//         OMEGAk = eph.OMEGA + (eph.OMEGAdot - omega_e) * tk - omega_e * eph.toe;
//     } else {
//         // GPS使用WGS84
//         OMEGAk = eph.OMEGA + (eph.OMEGAdot - omega_e) * tk - omega_e * eph.toe;
//     }

//     // 计算ECEF坐标
//     pos.X = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
//     pos.Y = xk * sin(OMEGAk) + yk * cos(ik) * cos(OMEGAk);
//     pos.Z = yk * sin(ik);
    
//     return pos;
// }

// 查找精密星历点
bool FindPrecisePoint(int prn, char satSystem, double time, int hour, int minute,
                     const vector<PreciseEphemerisPoint>& precisePoints,
                     PreciseEphemerisPoint& precisePoint) {
    for (const auto& point : precisePoints) {
        // 首先检查时间是否匹配
        if (point.hour == hour && point.minute == minute) {
            if (point.PRN == prn && point.satSystem == satSystem) {
                precisePoint = point;
                return true;
            }
        }
    }
    
    // 如果精确时间不匹配，查找最接近的
    double minTimeDiff = 1e9;
    bool found = false;
    
    for (const auto& point : precisePoints) {
        if (point.PRN == prn && point.satSystem == satSystem) {
            double timeDiff = fabs(point.time - time);
            if (timeDiff < minTimeDiff && timeDiff < 300.0) { // 限制在5分钟内
                minTimeDiff = timeDiff;
                precisePoint = point;
                found = true;
            }
        }
    }
    
    return found;
}

// GPS/BDS卫星位置计算（修正版）
SatellitePosition CalculateSatellitePosition(const EPHEMERISBLOCK& eph, double transmitTime, int hour, int minute) {
    SatellitePosition pos;
    pos.PRN = eph.PRN;
    pos.satSystem = eph.satSystem;
    pos.time = transmitTime;
    pos.hour = hour;
    pos.minute = minute;
    pos.tow = fmod(transmitTime, 86400.0);
    
    bool isBDS = (eph.satSystem == 'C');
    
    // 使用正确的常数
    double mu = isBDS ? MU_BDS : MU_GPS;
    double omega_e = isBDS ? OMEGA_E_BDS : OMEGA_E_GPS;
    
    // 关键修改：BDS时间系统处理
    // BDS使用BDT时间，与GPS时间有14秒偏差
    double timeOffset = isBDS ? 14.0 : 0.0;
    double adjustedTime = transmitTime - timeOffset;
    
    // 计算时间差 tk = t - toe（使用调整后的时间）
    double tk = adjustedTime - eph.toe;
    
    // 处理周跳（根据图片说明）
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;
    
    // 计算平均角速度
    double A = eph.SqrtA * eph.SqrtA;
    double n0 = sqrt(mu / (A * A * A));
    double n = n0 + eph.Deltan;
    
    // 计算平近点角
    double Mk = eph.M0 + n * tk;
    
    // 迭代计算偏近点角（开普勒方程）
    double Ek = Mk;
    for (int i = 0; i < 10; i++) {
        double Ek_prev = Ek;
        Ek = Mk + eph.e * sin(Ek);
        if (fabs(Ek - Ek_prev) < 1e-12) break;
    }
    
    // 计算真近点角
    double sin_vk = (sqrt(1 - eph.e * eph.e) * sin(Ek)) / (1 - eph.e * cos(Ek));
    double cos_vk = (cos(Ek) - eph.e) / (1 - eph.e * cos(Ek));
    double vk = atan2(sin_vk, cos_vk);
    
    // 计算纬度幅角
    double phik = vk + eph.omega;
    
    // 计算周期改正项
    double delta_uk = eph.Cus * sin(2 * phik) + eph.Cuc * cos(2 * phik);
    double delta_rk = eph.Crs * sin(2 * phik) + eph.Crc * cos(2 * phik);
    double delta_ik = eph.Cis * sin(2 * phik) + eph.Cic * cos(2 * phik);
    
    // 计算改正后的参数
    double uk = phik + delta_uk;
    double rk = A * (1 - eph.e * cos(Ek)) + delta_rk;
    double ik = eph.i0 + eph.IDOT * tk + delta_ik;
    
    // 计算轨道平面坐标
    double xk = rk * cos(uk);
    double yk = rk * sin(uk);
    
    // 关键修改：BDS卫星分类处理
    if (isBDS) {
        // 判断卫星类型：GEO卫星(PRN 1-5), MEO/IGSO卫星(PRN 6-35)
        bool isGEO = (eph.PRN >= 1 && eph.PRN <= 5);
        
        if (isGEO) {
            // GEO卫星特殊处理（根据图片中的GEO公式）
            // 计算历元升交点经度（惯性系）
            double OMEGAk_inertial = eph.OMEGA + eph.OMEGAdot * tk - omega_e * eph.toe;
            
            // 在惯性系中的坐标
            double X_GK = xk * cos(OMEGAk_inertial) - yk * cos(ik) * sin(OMEGAk_inertial);
            double Y_GK = xk * sin(OMEGAk_inertial) + yk * cos(ik) * cos(OMEGAk_inertial);
            double Z_GK = yk * sin(ik);
            
            // 坐标转换到BDCS系：绕Z轴旋转(ωe*tk)，再绕X轴旋转(-5°)
            double rotationAngle = omega_e * tk;
            double phi = -5.0 * PI / 180.0; // -5度转弧度
            
            // 绕Z轴旋转
            double X_temp = X_GK * cos(rotationAngle) + Y_GK * sin(rotationAngle);
            double Y_temp = -X_GK * sin(rotationAngle) + Y_GK * cos(rotationAngle);
            double Z_temp = Z_GK;
            
            // 绕X轴旋转(-5°)
            pos.Y = Y_temp * cos(phi) + Z_temp * sin(phi);
            pos.Z = -Y_temp * sin(phi) + Z_temp * cos(phi);
            pos.X = X_temp;
            
        } else {
            // MEO/IGSO卫星处理（与GPS类似但使用BDCS参数）
            double OMEGAk = eph.OMEGA + (eph.OMEGAdot - omega_e) * tk - omega_e * eph.toe;
            
            pos.X = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
            pos.Y = xk * sin(OMEGAk) + yk * cos(ik) * cos(OMEGAk);
            pos.Z = yk * sin(ik);
        }
    } else {
        // GPS卫星计算（保持不变）
        double OMEGAk = eph.OMEGA + (eph.OMEGAdot - omega_e) * tk - omega_e * eph.toe;
        
        pos.X = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
        pos.Y = xk * sin(OMEGAk) + yk * cos(ik) * cos(OMEGAk);
        pos.Z = yk * sin(ik);
    }
    
    return pos;
}

// 添加旋转矩阵计算函数
void RotateZ(double angle, double& x, double& y, double& z) {
    double cosA = cos(angle);
    double sinA = sin(angle);
    double x_temp = x * cosA + y * sinA;
    double y_temp = -x * sinA + y * cosA;
    x = x_temp;
    y = y_temp;
    // z保持不变
}

void RotateX(double angle, double& x, double& y, double& z) {
    double cosA = cos(angle);
    double sinA = sin(angle);
    double y_temp = y * cosA + z * sinA;
    double z_temp = -y * sinA + z * cosA;
    y = y_temp;
    z = z_temp;
    // x保持不变
}


// 比较广播星历与精密星历
vector<PositionComparison> CompareEphemeris(const vector<SatellitePosition>& broadcastPositions,
                                          const vector<PreciseEphemerisPoint>& precisePoints) {
    vector<PositionComparison> comparisons;
    
    for (const auto& bpos : broadcastPositions) {
        PreciseEphemerisPoint ppos;
        
        if (FindPrecisePoint(bpos.PRN, bpos.satSystem, bpos.time, bpos.hour, bpos.minute, precisePoints, ppos)) {
            PositionComparison comp;
            comp.PRN = bpos.PRN;
            comp.satSystem = bpos.satSystem;
            comp.time = bpos.time;
            comp.hour = bpos.hour;
            comp.minute = bpos.minute;
            
            comp.broadcastX = bpos.X;
            comp.broadcastY = bpos.Y;
            comp.broadcastZ = bpos.Z;
            
            comp.preciseX = ppos.X;
            comp.preciseY = ppos.Y;
            comp.preciseZ = ppos.Z;
            
            comp.dX = bpos.X - ppos.X;
            comp.dY = bpos.Y - ppos.Y;
            comp.dZ = bpos.Z - ppos.Z;
            
            comp.distance = sqrt(comp.dX * comp.dX + comp.dY * comp.dY + comp.dZ * comp.dZ);
            
            comparisons.push_back(comp);
        }
    }
    
    return comparisons;
}

// 保存比较结果到TXT文件
void SaveComparisonResults(const vector<PositionComparison>& comparisons, const char* filename) {
    ofstream outfile(filename);
    
    if (!outfile.is_open()) {
        cerr << "Error: Cannot open output file " << filename << endl;
        return;
    }
    
    // 写入文件头
    outfile << "Broadcast vs Precise Ephemeris Comparison Results (GPS and BDS only)\n";
    outfile << "======================================================================\n\n";
    outfile << "PRN  Sys  Time(HH:MM)   Broadcast X(m)       Broadcast Y(m)       Broadcast Z(m)       "
            << "Precise X(m)         Precise Y(m)         Precise Z(m)         "
            << "dX(m)               dY(m)               dZ(m)               Distance(m)\n";
    outfile << string(200, '-') << "\n";
    
    // 设置输出格式
    outfile.precision(3);
    outfile << fixed;
    
    // 按卫星系统、时间和PRN排序
    vector<PositionComparison> sortedComparisons = comparisons;
    sort(sortedComparisons.begin(), sortedComparisons.end(), 
         [](const PositionComparison& a, const PositionComparison& b) {
             if (a.satSystem != b.satSystem) return a.satSystem < b.satSystem;  // GPS在前
             if (a.time == b.time) return a.PRN < b.PRN;
             return a.time < b.time;
         });
    
    // 写入数据
    int gpsCount = 0, bdsCount = 0;
    for (const auto& comp : sortedComparisons) {
        char timeStr[10];
        sprintf(timeStr, "%02d:%02d", comp.hour, comp.minute);
        
        outfile << setw(3) << comp.PRN << "   " << comp.satSystem << "   " 
                << timeStr << "   "
                << setw(15) << comp.broadcastX << "   " << setw(15) << comp.broadcastY << "   " 
                << setw(15) << comp.broadcastZ << "   "
                << setw(15) << comp.preciseX << "   " << setw(15) << comp.preciseY << "   " 
                << setw(15) << comp.preciseZ << "   "
                << setw(15) << comp.dX << "   " << setw(15) << comp.dY << "   " 
                << setw(15) << comp.dZ << "   "
                << setw(15) << comp.distance << "\n";
        
        if (comp.satSystem == 'G') gpsCount++;
        else if (comp.satSystem == 'C') bdsCount++;
    }
    
    // 计算统计信息
    if (!sortedComparisons.empty()) {
        double sumDistance = 0;
        double maxDistance = 0;
        double minDistance = 1e9;
        
        map<int, vector<double>> gpsDistances;  // GPS卫星统计
        map<int, vector<double>> bdsDistances;  // BDS卫星统计
        
        for (const auto& comp : sortedComparisons) {
            sumDistance += comp.distance;
            if (comp.distance > maxDistance) maxDistance = comp.distance;
            if (comp.distance < minDistance) minDistance = comp.distance;
            
            if (comp.satSystem == 'G') {
                gpsDistances[comp.PRN].push_back(comp.distance);
            } else if (comp.satSystem == 'C') {
                bdsDistances[comp.PRN].push_back(comp.distance);
            }
        }
        
        double avgDistance = sumDistance / sortedComparisons.size();
        
        outfile << "\n\nOverall Statistics:\n";
        outfile << "===================================================\n";
        outfile << "Total comparisons: " << sortedComparisons.size() << " (GPS: " << gpsCount << ", BDS: " << bdsCount << ")\n";
        outfile << "Average position error: " << avgDistance << " m\n";
        outfile << "Maximum position error: " << maxDistance << " m\n";
        outfile << "Minimum position error: " << minDistance << " m\n";
        
        // GPS卫星统计
        if (!gpsDistances.empty()) {
            outfile << "\n\nGPS Satellite Statistics:\n";
            outfile << "PRN  Count  Avg.Error(m)  Max.Error(m)  Min.Error(m)\n";
            outfile << "---------------------------------------------------\n";
            
            for (const auto& entry : gpsDistances) {
                int prn = entry.first;
                const vector<double>& dists = entry.second;
                
                double satSum = 0, satMax = 0, satMin = 1e9;
                for (double d : dists) {
                    satSum += d;
                    if (d > satMax) satMax = d;
                    if (d < satMin) satMin = d;
                }
                double satAvg = satSum / dists.size();
                
                outfile << setw(3) << prn << "   " << setw(5) << dists.size() << "   "
                        << setw(12) << satAvg << "   " << setw(12) << satMax << "   "
                        << setw(12) << satMin << "\n";
            }
        }
        
        // BDS卫星统计
        if (!bdsDistances.empty()) {
            outfile << "\n\nBDS Satellite Statistics:\n";
            outfile << "PRN  Count  Avg.Error(m)  Max.Error(m)  Min.Error(m)\n";
            outfile << "---------------------------------------------------\n";
            
            for (const auto& entry : bdsDistances) {
                int prn = entry.first;
                const vector<double>& dists = entry.second;
                
                double satSum = 0, satMax = 0, satMin = 1e9;
                for (double d : dists) {
                    satSum += d;
                    if (d > satMax) satMax = d;
                    if (d < satMin) satMin = d;
                }
                double satAvg = satSum / dists.size();
                
                outfile << setw(3) << prn << "   " << setw(5) << dists.size() << "   "
                        << setw(12) << satAvg << "   " << setw(12) << satMax << "   "
                        << setw(12) << satMin << "\n";
            }
        }
    }
    
    outfile.close();
    cout << "Comparison results saved to " << filename << endl;
}

// 生成MATLAB绘图代码（只绘制GPS卫星）
void GenerateMatlabPlotCode(const vector<SatellitePosition>& positions, const char* filename) {
    ofstream matlabFile(filename);
    
    matlabFile << "% MATLAB code for GPS satellite trajectory plotting\n";
    matlabFile << "% Generated from C++ program\n";
    matlabFile << "clear; close all; clc;\n\n";
    
    // 按卫星PRN组织数据（只处理GPS）
    map<int, vector<double>> gpsX, gpsY, gpsZ;
    map<int, vector<double>> bdsX, bdsY, bdsZ;
    map<int, vector<double>> gpsTime;
    
    for (const auto& pos : positions) {
        if (pos.satSystem == 'G') {
            gpsX[pos.PRN].push_back(pos.X);
            gpsY[pos.PRN].push_back(pos.Y);
            gpsZ[pos.PRN].push_back(pos.Z);
            gpsTime[pos.PRN].push_back(pos.tow / 3600.0); // 转换为小时
        } else if (pos.satSystem == 'C') {
            bdsX[pos.PRN].push_back(pos.X);
            bdsY[pos.PRN].push_back(pos.Y);
            bdsZ[pos.PRN].push_back(pos.Z);
        }
    }
    
    // 生成卫星数据
    matlabFile << "% GPS satellite positions data\n";
    
    // GPS卫星数据
    matlabFile << "% GPS satellites\n";
    for (const auto& entry : gpsX) {
        int prn = entry.first;
        matlabFile << "gps" << prn << "_X = [";
        for (double x : gpsX[prn]) matlabFile << x << " ";
        matlabFile << "];\n";
        
        matlabFile << "gps" << prn << "_Y = [";
        for (double y : gpsY[prn]) matlabFile << y << " ";
        matlabFile << "];\n";
        
        matlabFile << "gps" << prn << "_Z = [";
        for (double z : gpsZ[prn]) matlabFile << z << " ";
        matlabFile << "];\n\n";
    }
    
    // BDS卫星数据（可选）
    if (!bdsX.empty()) {
        matlabFile << "% BDS satellites (optional)\n";
        for (const auto& entry : bdsX) {
            int prn = entry.first;
            matlabFile << "bds" << prn << "_X = [";
            for (double x : bdsX[prn]) matlabFile << x << " ";
            matlabFile << "];\n";
            
            matlabFile << "bds" << prn << "_Y = [";
            for (double y : bdsY[prn]) matlabFile << y << " ";
            matlabFile << "];\n";
            
            matlabFile << "bds" << prn << "_Z = [";
            for (double z : bdsZ[prn]) matlabFile << z << " ";
            matlabFile << "];\n\n";
        }
    }
    
    // 绘图代码 - 主要绘制GPS卫星
    matlabFile << "% Plot GPS satellite trajectories\n";
    matlabFile << "figure('Position', [100, 100, 1400, 900]);\n";
    matlabFile << "hold on; grid on; box on;\n";
    
    // 设置坐标轴
    matlabFile << "xlabel('X (m)');\n";
    matlabFile << "ylabel('Y (m)');\n";
    matlabFile << "zlabel('Z (m)');\n";
    matlabFile << "title('GPS Satellite Trajectories in ECEF Coordinate System', 'FontSize', 14);\n";
    
    // 绘制地球
    matlabFile << "% Draw Earth sphere\n";
    matlabFile << "RE = 6378137; % Earth radius in meters\n";
    matlabFile << "[Xe, Ye, Ze] = sphere(50);\n";
    matlabFile << "surf(Xe*RE, Ye*RE, Ze*RE, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);\n";
    
    // 绘制GPS卫星轨道
    if (!gpsX.empty()) {
        matlabFile << "\n% Plot GPS satellite trajectories\n";
        matlabFile << "gpsColors = lines(" << gpsX.size() << ");\n";
        
        int colorIdx = 0;
        for (const auto& entry : gpsX) {
            int prn = entry.first;
            matlabFile << "plot3(gps" << prn << "_X, gps" << prn << "_Y, gps" << prn << "_Z, ...\n";
            matlabFile << "    'Color', gpsColors(" << (colorIdx + 1) << ",:), 'LineWidth', 2.0, ...\n";
            matlabFile << "    'DisplayName', sprintf('G%02d', " << prn << "));\n";
            
            // 标记起始点
            if (!gpsX[prn].empty()) {
                matlabFile << "plot3(gps" << prn << "_X(1), gps" << prn << "_Y(1), gps" << prn << "_Z(1), ...\n";
                matlabFile << "    'o', 'MarkerSize', 8, 'MarkerFaceColor', gpsColors(" << (colorIdx + 1) << ",:), ...\n";
                matlabFile << "    'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');\n";
            }
            
            colorIdx++;
        }
    }
    
    // 可选：绘制BDS卫星轨道（用虚线）
    if (!bdsX.empty()) {
        matlabFile << "\n% Plot BDS satellite trajectories (optional, dashed lines)\n";
        matlabFile << "bdsColor = [0.5 0.5 0.5]; % Gray color for BDS\n";
        for (const auto& entry : bdsX) {
            int prn = entry.first;
            matlabFile << "plot3(bds" << prn << "_X, bds" << prn << "_Y, bds" << prn << "_Z, ...\n";
            matlabFile << "    '--', 'Color', bdsColor, 'LineWidth', 1.5, ...\n";
            matlabFile << "    'DisplayName', sprintf('C%02d', " << prn << "));\n";
        }
    }
    
    // 设置图形属性
    matlabFile << "\n% Set plot properties\n";
    matlabFile << "axis equal;\n";
    matlabFile << "view(45, 30);\n";
    matlabFile << "legend('Location', 'best', 'FontSize', 10);\n";
    matlabFile << "set(gca, 'FontSize', 12);\n";
    matlabFile << "grid on;\n";
    
    // 添加光照效果
    matlabFile << "light('Position', [1 1 1], 'Style', 'infinite');\n";
    matlabFile << "lighting gouraud;\n";
    
    // 保存图形
    matlabFile << "\n% Save figure\n";
    matlabFile << "saveas(gcf, 'gps_satellite_trajectories_3d.png');\n";
    matlabFile << "disp('GPS 3D trajectory plot saved as gps_satellite_trajectories_3d.png');\n";
    
    // 创建误差分析图
    matlabFile << "\n% Create error analysis plot\n";
    matlabFile << "figure('Position', [200, 200, 1200, 500]);\n";
    
    // 读取比较结果文件
    matlabFile << "% Read comparison results\n";
    matlabFile << "fid = fopen('ephemeris_comparison_results.txt', 'r');\n";
    matlabFile << "if fid ~= -1\n";
    matlabFile << "    fclose(fid);\n";
    matlabFile << "    % You can add code here to analyze the error data\n";
    matlabFile << "    disp('Comparison results available for analysis');\n";
    matlabFile << "end\n";
    
    matlabFile.close();
    cout << "MATLAB plotting code generated: " << filename << endl;
}

int main() {
    // 设置文件路径
    const char* broadcastFile = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/brdm3350.19p";
    const char* preciseFile = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/WUM0MGXFIN_20193350000_01D_15M_ORB.SP3";
    
    vector<EPHEMERISBLOCK> ephemerisBlocks;
    
    cout << "Reading broadcast ephemeris file..." << endl;
    if (ReadBrodcastEphemeris(broadcastFile, ephemerisBlocks) != 0) {
        cerr << "Failed to read broadcast ephemeris file!" << endl;
        return -1;
    }
    
    if (ephemerisBlocks.empty()) {
        cerr << "No GPS or BDS ephemeris data found!" << endl;
        return -1;
    }
    
    // 读取精密星历文件
    cout << "Reading precise ephemeris file..." << endl;
    vector<PreciseEphemerisPoint> precisePoints = ReadPreciseEphemeris(preciseFile);
    
    if (precisePoints.empty()) {
        cerr << "No precise ephemeris data found!" << endl;
        return -1;
    }
    
    // 收集所有GPS和BDS卫星PRN
    map<pair<int, char>, bool> allSatellites;
    for (const auto& point : precisePoints) {
        if (point.satSystem == 'G' || point.satSystem == 'C') {
            allSatellites[{point.PRN, point.satSystem}] = true;
        }
    }
    
    cout << "\nFound " << allSatellites.size() << " GPS/BDS satellites in precise ephemeris." << endl;
    
    // 计算卫星位置（只计算精密星历中存在的时间点）
    vector<SatellitePosition> broadcastPositions;
    
    // 从精密星历中提取所有时间点
    map<pair<int, int>, bool> timePoints;  // hour, minute
    for (const auto& point : precisePoints) {
        if (point.satSystem == 'G' || point.satSystem == 'C') {
            timePoints[{point.hour, point.minute}] = true;
        }
    }
    
    cout << "\nCalculating broadcast ephemeris positions..." << endl;
    
    int calculatedCount = 0;
    for (const auto& timePoint : timePoints) {
        int hour = timePoint.first.first;
        int minute = timePoint.first.second;
        double transmitTime = hour * 3600.0 + minute * 60.0;
        
        // 为每个GPS/BDS卫星计算位置
        for (const auto& satEntry : allSatellites) {
            int prn = satEntry.first.first;
            char satSystem = satEntry.first.second;
            
            const EPHEMERISBLOCK* eph = FindClosestEphemeris(prn, satSystem, transmitTime, ephemerisBlocks);
            if (eph) {
                SatellitePosition pos = CalculateSatellitePosition(*eph, transmitTime, hour, minute);
                broadcastPositions.push_back(pos);
                calculatedCount++;
            }
        }
    }
    
    cout << "Calculated " << calculatedCount << " broadcast ephemeris positions." << endl;
    
    // 输出部分结果
    if (!broadcastPositions.empty()) {
        cout << "\nFirst 10 calculated positions (GPS/BDS only):" << endl;
        cout << "PRN Sys Time(HH:MM)    X(m)            Y(m)            Z(m)" << endl;
        cout.precision(3);
        cout << fixed;
        
        int count = 0;
        for (const auto& pos : broadcastPositions) {
            if (count >= 10) break;
            cout << setw(3) << pos.PRN << "   " << pos.satSystem << "   "
                 << pos.hour << ":" << setw(2) << setfill('0') << pos.minute << setfill(' ') << "   "
                 << setw(12) << pos.X << "   " << setw(12) << pos.Y << "   " << setw(12) << pos.Z << endl;
            count++;
        }
    }
    
    // 比较广播星历与精密星历
    cout << "\nComparing broadcast and precise ephemeris..." << endl;
    vector<PositionComparison> comparisons = CompareEphemeris(broadcastPositions, precisePoints);
    
    if (!comparisons.empty()) {
        cout << "Found " << comparisons.size() << " matching points for comparison." << endl;
        
        // 保存比较结果
        SaveComparisonResults(comparisons, "ephemeris_comparison_results.txt");
        
        // 生成MATLAB绘图代码
        GenerateMatlabPlotCode(broadcastPositions, "gps_trajectories_3d.m");
        cout << "\nMATLAB plotting code generated: gps_trajectories_3d.m" << endl;
        cout << "Run this file in MATLAB to generate the 3D trajectory plot." << endl;
    } else {
        cout << "Warning: No matching points found for comparison!" << endl;
        cout << "Make sure the broadcast and precise ephemeris cover the same time period and satellites." << endl;
    }
    
    cout << "\nProgram completed successfully!" << endl;
    return 0;
}