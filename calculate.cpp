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
    
    // 钟差参数
    double a0, a1, a2;
    
    // 轨道参数
    double IODE, Crs, Deltan, M0;        // ORBIT - 1
    double Cuc, e, Cus, SqrtA;           // ORBIT - 2  
    double Toe, Cic, OMEGA, Cis;         // ORBIT - 3
    double i0, Crc, omega, OMEGAdot;     // ORBIT - 4
    double IDOT, GpsWeekNumber, L2C, L2P; // ORBIT - 5
    double SatAccuracy, SatHealth, TGD, IODC; // ORBIT - 6
    
    // BDS特有参数
    double AODE;  // BDS星历数据龄期
    double AODC;  // BDS钟差数据龄期
};

struct SatellitePosition {
    int PRN;
    char satSystem;
    double X, Y, Z;  // ECEF坐标 (m)
    double time;      // 时间 (GPS秒)
    int hour, minute; // 计算时间的小时和分钟
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

// 读取广播星历文件
int ReadBrodcastEphemeris(const char* filename, int &EphemerisBlockNum, EPHEMERISBLOCK* &pGpsEphemeris) {
    FILE* pfEph = fopen(filename, "r");
    if (!pfEph) {
        cerr << "Error: Cannot open file " << filename << endl;
        return -1;
    }
    
    // 第一遍：计算星历块数量
    char strLine[256];
    int headLineNum = 0;
    
    // 跳过文件头
    while (fgets(strLine, sizeof(strLine), pfEph)) {
        headLineNum++;
        if (strstr(strLine, "END OF HEADER")) {
            break;
        }
    }
    
    // 计算星历块数量
    EphemerisBlockNum = 0;
    while (fgets(strLine, sizeof(strLine), pfEph)) {
        if (strLine[0] == 'G' || strLine[0] == 'C' || strLine[0] == 'R' || 
            strLine[0] == 'E' || strLine[0] == 'J') {
            EphemerisBlockNum++;
            // 跳过该卫星的后续行
            for (int i = 0; i < 7; i++) {
                fgets(strLine, sizeof(strLine), pfEph);
            }
        }
    }
    
    // 分配内存
    pGpsEphemeris = new EPHEMERISBLOCK[EphemerisBlockNum];
    
    // 第二遍：读取数据
    fseek(pfEph, 0, SEEK_SET);
    
    // 跳过文件头
    while (fgets(strLine, sizeof(strLine), pfEph)) {
        if (strstr(strLine, "END OF HEADER")) {
            break;
        }
    }
    
    int blockIndex = 0;
    while (fgets(strLine, sizeof(strLine), pfEph) && blockIndex < EphemerisBlockNum) {
        char satType;
        int prn, year, month, day, hour, minute;
        double second, a0, a1, a2;
        
        if (sscanf(strLine, "%c%2d %d %d %d %d %d %lf %lf %lf %lf", 
                   &satType, &prn, &year, &month, &day, &hour, &minute, &second, &a0, &a1, &a2) >= 9) {
            
            pGpsEphemeris[blockIndex].satSystem = satType;
            pGpsEphemeris[blockIndex].PRN = prn;
            pGpsEphemeris[blockIndex].year = year;
            pGpsEphemeris[blockIndex].month = month;
            pGpsEphemeris[blockIndex].day = day;
            pGpsEphemeris[blockIndex].hour = hour;
            pGpsEphemeris[blockIndex].minute = minute;
            pGpsEphemeris[blockIndex].second = second;
            pGpsEphemeris[blockIndex].a0 = a0;
            pGpsEphemeris[blockIndex].a1 = a1;
            pGpsEphemeris[blockIndex].a2 = a2;
            
            // 读取轨道参数
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &pGpsEphemeris[blockIndex].IODE, 
                   &pGpsEphemeris[blockIndex].Crs, 
                   &pGpsEphemeris[blockIndex].Deltan, 
                   &pGpsEphemeris[blockIndex].M0);
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &pGpsEphemeris[blockIndex].Cuc, 
                   &pGpsEphemeris[blockIndex].e, 
                   &pGpsEphemeris[blockIndex].Cus, 
                   &pGpsEphemeris[blockIndex].SqrtA);
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &pGpsEphemeris[blockIndex].Toe, 
                   &pGpsEphemeris[blockIndex].Cic, 
                   &pGpsEphemeris[blockIndex].OMEGA, 
                   &pGpsEphemeris[blockIndex].Cis);
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &pGpsEphemeris[blockIndex].i0, 
                   &pGpsEphemeris[blockIndex].Crc, 
                   &pGpsEphemeris[blockIndex].omega, 
                   &pGpsEphemeris[blockIndex].OMEGAdot);
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &pGpsEphemeris[blockIndex].IDOT, 
                   &pGpsEphemeris[blockIndex].L2C, 
                   &pGpsEphemeris[blockIndex].GpsWeekNumber, 
                   &pGpsEphemeris[blockIndex].L2P);
            
            fgets(strLine, sizeof(strLine), pfEph);
            sscanf(strLine, "%lf %lf %lf %lf", 
                   &pGpsEphemeris[blockIndex].SatAccuracy, 
                   &pGpsEphemeris[blockIndex].SatHealth, 
                   &pGpsEphemeris[blockIndex].TGD, 
                   &pGpsEphemeris[blockIndex].IODC);
            
            // 跳过最后一行
            fgets(strLine, sizeof(strLine), pfEph);
            
            blockIndex++;
        }
    }
    
    fclose(pfEph);
    return 0;
}

// 读取精密星历文件 (SP3格式)
vector<PreciseEphemerisPoint> ReadPreciseEphemeris(const char* filename) {
    vector<PreciseEphemerisPoint> precisePoints;
    
    ifstream sp3file(filename);
    if (!sp3file.is_open()) {
        cerr << "Error: Cannot open precise ephemeris file " << filename << endl;
        return precisePoints;
    }
    
    string line;
    int lineCount = 0;
    int epochCount = 0;
    int year, month, day, hour, minute;
    double second;
    double gpsWeekSecond = 0;
    
    // 读取文件头
    while (getline(sp3file, line)) {
        lineCount++;
        
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
                epochCount++;
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
            PreciseEphemerisPoint point;
            point.PRN = prn;
            point.satSystem = satType;
            point.X = x * 1000.0;  // km转换为m
            point.Y = y * 1000.0;
            point.Z = z * 1000.0;
            point.time = gpsWeekSecond;
            point.hour = hour;
            point.minute = minute;
            
            precisePoints.push_back(point);
        }
    }
    
    sp3file.close();
    
    cout << "Read " << precisePoints.size() << " precise ephemeris points from " << filename << endl;
    return precisePoints;
}

// GPS卫星位置计算
SatellitePosition CalculateGpsSatellitePosition(const EPHEMERISBLOCK& eph, double transmitTime, int hour, int minute) {
    SatellitePosition pos;
    pos.PRN = eph.PRN;
    pos.satSystem = 'G';
    pos.time = transmitTime;
    pos.hour = hour;
    pos.minute = minute;
    
    double mu = MU_GPS;
    double omega_e = OMEGA_E_GPS;
    
    // 计算时间差
    double weekSecond;
    Calendar2GpsTime(eph.year, eph.month, eph.day, eph.hour, eph.minute, eph.second, weekSecond);
    double toe = eph.Toe;
    double tk = transmitTime - toe;
    
    // 处理周跳
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;
    
    // 计算平均角速度
    double A = eph.SqrtA * eph.SqrtA;
    double n0 = sqrt(mu / (A * A * A));
    double n = n0 + eph.Deltan;
    
    // 计算平近点角
    double Mk = eph.M0 + n * tk;
    
    // 迭代计算偏近点角
    double Ek = Mk;
    double Ek_prev;
    for (int i = 0; i < 10; i++) {
        Ek_prev = Ek;
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
    
    // 计算升交点经度
    double OMEGAk = eph.OMEGA + (eph.OMEGAdot - omega_e) * tk - omega_e * toe;
    
    // 计算ECEF坐标
    pos.X = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
    pos.Y = xk * sin(OMEGAk) + yk * cos(ik) * cos(OMEGAk);
    pos.Z = yk * sin(ik);
    
    return pos;
}

// BDS卫星位置计算
SatellitePosition CalculateBdsSatellitePosition(const EPHEMERISBLOCK& eph, double transmitTime, int hour, int minute) {
    SatellitePosition pos;
    pos.PRN = eph.PRN;
    pos.satSystem = 'C';
    pos.time = transmitTime;
    pos.hour = hour;
    pos.minute = minute;
    
    double mu = MU_BDS;
    double omega_e = OMEGA_E_BDS;
    
    // 计算时间差
    double weekSecond;
    Calendar2GpsTime(eph.year, eph.month, eph.day, eph.hour, eph.minute, eph.second, weekSecond);
    double toe = eph.Toe;
    double tk = transmitTime - toe;
    
    // 处理周跳
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;
    
    // 计算平均角速度
    double A = eph.SqrtA * eph.SqrtA;
    double n0 = sqrt(mu / (A * A * A));
    double n = n0 + eph.Deltan;
    
    // 计算平近点角
    double Mk = eph.M0 + n * tk;
    
    // 迭代计算偏近点角
    double Ek = Mk;
    double Ek_prev;
    for (int i = 0; i < 10; i++) {
        Ek_prev = Ek;
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
    
    // 计算升交点经度
    double OMEGAk = eph.OMEGA + (eph.OMEGAdot - omega_e) * tk - omega_e * toe;
    
    // 计算ECEF坐标
    pos.X = xk * cos(OMEGAk) - yk * cos(ik) * sin(OMEGAk);
    pos.Y = xk * sin(OMEGAk) + yk * cos(ik) * cos(OMEGAk);
    pos.Z = yk * sin(ik);
    
    return pos;
}

// 查找最接近的精密星历点
bool FindClosestPrecisePoint(int prn, char satSystem, double time, 
                           const vector<PreciseEphemerisPoint>& precisePoints,
                           PreciseEphemerisPoint& closestPoint) {
    double minTimeDiff = 1e9;
    bool found = false;
    
    for (const auto& point : precisePoints) {
        if (point.PRN == prn && point.satSystem == satSystem) {
            double timeDiff = fabs(point.time - time);
            if (timeDiff < minTimeDiff && timeDiff < 300.0) { // 限制在5分钟内
                minTimeDiff = timeDiff;
                closestPoint = point;
                found = true;
            }
        }
    }
    
    return found;
}

// 比较广播星历与精密星历
vector<PositionComparison> CompareEphemeris(const vector<SatellitePosition>& broadcastPositions,
                                          const vector<PreciseEphemerisPoint>& precisePoints) {
    vector<PositionComparison> comparisons;
    
    for (const auto& bpos : broadcastPositions) {
        PreciseEphemerisPoint ppos;
        
        if (FindClosestPrecisePoint(bpos.PRN, bpos.satSystem, bpos.time, precisePoints, ppos)) {
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
    outfile << "Broadcast vs Precise Ephemeris Comparison Results\n";
    outfile << "===================================================\n\n";
    outfile << "PRN Sys  Time(HH:MM)   Broadcast X(m)       Broadcast Y(m)       Broadcast Z(m)       "
            << "Precise X(m)         Precise Y(m)         Precise Z(m)         "
            << "dX(m)               dY(m)               dZ(m)               Distance(m)\n";
    outfile << string(200, '-') << "\n";
    
    // 设置输出格式
    outfile.precision(3);
    outfile << fixed;
    
    // 写入数据
    for (const auto& comp : comparisons) {
        char timeStr[10];
        sprintf(timeStr, "%02d:%02d", comp.hour, comp.minute);
        
        outfile << comp.PRN << "   " << comp.satSystem << "   " 
                << timeStr << "   "
                << comp.broadcastX << "   " << comp.broadcastY << "   " << comp.broadcastZ << "   "
                << comp.preciseX << "   " << comp.preciseY << "   " << comp.preciseZ << "   "
                << comp.dX << "   " << comp.dY << "   " << comp.dZ << "   "
                << comp.distance << "\n";
    }
    
    // 计算统计信息
    if (!comparisons.empty()) {
        double sumDistance = 0;
        double maxDistance = 0;
        double minDistance = 1e9;
        
        for (const auto& comp : comparisons) {
            sumDistance += comp.distance;
            if (comp.distance > maxDistance) maxDistance = comp.distance;
            if (comp.distance < minDistance) minDistance = comp.distance;
        }
        
        double avgDistance = sumDistance / comparisons.size();
        
        outfile << "\n\nStatistics:\n";
        outfile << "Total comparisons: " << comparisons.size() << "\n";
        outfile << "Average distance error: " << avgDistance << " m\n";
        outfile << "Maximum distance error: " << maxDistance << " m\n";
        outfile << "Minimum distance error: " << minDistance << " m\n";
    }
    
    outfile.close();
    cout << "Comparison results saved to " << filename << endl;
}

// 生成MATLAB绘图代码 - 改进版，每颗卫星单独轨道
void GenerateMatlabPlotCode(const vector<SatellitePosition>& positions, const char* filename) {
    ofstream matlabFile(filename);
    
    matlabFile << "% MATLAB code for satellite trajectory plotting\n";
    matlabFile << "% Generated from C++ program\n";
    matlabFile << "clear; close all; clc;\n\n";
    
    // 按卫星PRN组织数据
    map<int, vector<double>> gpsX, gpsY, gpsZ;
    map<int, vector<double>> bdsX, bdsY, bdsZ;
    
    for (const auto& pos : positions) {
        if (pos.satSystem == 'G') {
            gpsX[pos.PRN].push_back(pos.X);
            gpsY[pos.PRN].push_back(pos.Y);
            gpsZ[pos.PRN].push_back(pos.Z);
        } else if (pos.satSystem == 'C') {
            bdsX[pos.PRN].push_back(pos.X);
            bdsY[pos.PRN].push_back(pos.Y);
            bdsZ[pos.PRN].push_back(pos.Z);
        }
    }
    
    // 生成卫星数据
    matlabFile << "% Satellite positions data\n";
    
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
    
    // BDS卫星数据
    matlabFile << "% BDS satellites\n";
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
    
    // 绘图代码 - 类似附件的风格
    matlabFile << "% Plot satellite trajectories\n";
    matlabFile << "figure('Position', [100, 100, 1400, 900]);\n";
    matlabFile << "hold on; grid on; box on;\n";
    
    // 设置坐标轴范围和标签
    matlabFile << "xlabel('X (m)');\n";
    matlabFile << "ylabel('Y (m)');\n";
    matlabFile << "zlabel('Z (m)');\n";
    matlabFile << "title('GPS Satellite Trajectories in ECEF Coordinate System', 'FontSize', 14);\n";
    
    // 绘制地球
    matlabFile << "% Draw Earth sphere\n";
    matlabFile << "RE = 6378137; % Earth radius in meters\n";
    matlabFile << "[Xe, Ye, Ze] = sphere(50);\n";
    matlabFile << "surf(Xe*RE, Ye*RE, Ze*RE, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);\n";
    
    // 使用不同的颜色和线型绘制每颗GPS卫星轨道
    matlabFile << "% Plot GPS satellite trajectories\n";
    matlabFile << "colors = lines(" << gpsX.size() << ");\n";
    
    int colorIdx = 0;
    for (const auto& entry : gpsX) {
        int prn = entry.first;
        matlabFile << "plot3(gps" << prn << "_X, gps" << prn << "_Y, gps" << prn << "_Z, ...\n";
        matlabFile << "    'Color', colors(" << (colorIdx + 1) << ",:), 'LineWidth', 2.0, 'DisplayName', 'G" << prn << "');\n";
        
        // 标记起始点
        if (!gpsX[prn].empty()) {
            matlabFile << "plot3(gps" << prn << "_X(1), gps" << prn << "_Y(1), gps" << prn << "_Z(1), ...\n";
            matlabFile << "    'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(" << (colorIdx + 1) << ",:), ...\n";
            matlabFile << "    'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');\n";
        }
        
        colorIdx++;
    }
    
    // 如果BDS卫星存在，用不同样式绘制
    if (!bdsX.empty()) {
        matlabFile << "\n% Plot BDS satellite trajectories (dashed lines)\n";
        for (const auto& entry : bdsX) {
            int prn = entry.first;
            matlabFile << "plot3(bds" << prn << "_X, bds" << prn << "_Y, bds" << prn << "_Z, ...\n";
            matlabFile << "    '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'C" << prn << "');\n";
        }
    }
    
    // 设置图形属性
    matlabFile << "\n% Set plot properties\n";
    matlabFile << "axis equal;\n";
    matlabFile << "view(45, 30); % 设置视角\n";
    matlabFile << "legend('Location', 'best', 'FontSize', 10);\n";
    matlabFile << "set(gca, 'FontSize', 12);\n";
    
    // 添加网格和光照效果
    matlabFile << "grid on;\n";
    matlabFile << "light('Position', [1 1 1], 'Style', 'infinite');\n";
    matlabFile << "lighting gouraud;\n";
    matlabFile << "material dull;\n";
    
    // 保存图形
    matlabFile << "\n% Save figure\n";
    matlabFile << "saveas(gcf, 'satellite_trajectories_3d.png');\n";
    matlabFile << "disp('3D trajectory plot saved as satellite_trajectories_3d.png');\n";
    
    // 生成统计数据图
    matlabFile << "\n% Create statistics figure\n";
    matlabFile << "figure('Position', [200, 200, 1000, 600]);\n";
    matlabFile << "subplot(2,2,1);\n";
    
    // 卫星数量统计
    matlabFile << "satCounts = [" << gpsX.size() << " " << bdsX.size() << "];\n";
    matlabFile << "bar(satCounts);\n";
    matlabFile << "set(gca, 'XTickLabel', {'GPS', 'BDS'});\n";
    matlabFile << "ylabel('Number of Satellites');\n";
    matlabFile << "title('Satellite Count by System');\n";
    matlabFile << "grid on;\n";
    
    matlabFile.close();
    cout << "MATLAB plotting code generated: " << filename << endl;
}

int main() {
    // 读取广播星历文件
    const char* broadcastFile = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/brdm3350.19p";
    int ephemerisBlockNum = 0;
    EPHEMERISBLOCK* ephemerisBlocks = nullptr;
    
    cout << "Reading broadcast ephemeris file..." << endl;
    if (ReadBrodcastEphemeris(broadcastFile, ephemerisBlockNum, ephemerisBlocks) != 0) {
        cerr << "Failed to read broadcast ephemeris file!" << endl;
        return -1;
    }
    cout << "Successfully read " << ephemerisBlockNum << " ephemeris blocks." << endl;
    
    // 读取精密星历文件
    const char* preciseFile = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/WUM0MGXFIN_20193350000_01D_15M_ORB.SP3";
    vector<PreciseEphemerisPoint> precisePoints = ReadPreciseEphemeris(preciseFile);
    
    // 计算当天任意时刻的卫星位置
    vector<SatellitePosition> allPositions;
    int calculationDate[] = {2019, 12, 1};  // 年月日
    
    // 计算从00:00:00到23:59:59，每隔15分钟计算一次（与精密星历时间间隔匹配）
    cout << "\nCalculating satellite positions..." << endl;
    for (int hour = 0; hour < 24; hour++) {
        for (int minute = 0; minute < 60; minute += 15) {
            double transmitTime = hour * 3600.0 + minute * 60.0;
            
            for (int i = 0; i < ephemerisBlockNum; i++) {
                SatellitePosition pos;
                
                if (ephemerisBlocks[i].satSystem == 'G') {
                    pos = CalculateGpsSatellitePosition(ephemerisBlocks[i], transmitTime, hour, minute);
                } else if (ephemerisBlocks[i].satSystem == 'C') {
                    pos = CalculateBdsSatellitePosition(ephemerisBlocks[i], transmitTime, hour, minute);
                } else {
                    continue;  // 跳过其他卫星系统
                }
                
                allPositions.push_back(pos);
            }
        }
    }
    
    // 输出部分结果
    cout << "\nCalculated positions for " << allPositions.size() << " satellite-time points." << endl;
    cout << "First 10 positions:" << endl;
    cout << "PRN Sys Time(HH:MM)    X(m)            Y(m)            Z(m)" << endl;
    cout.precision(3);
    cout << fixed;
    
    for (int i = 0; i < min(10, (int)allPositions.size()); i++) {
        const auto& pos = allPositions[i];
        cout << pos.PRN << "   " << pos.satSystem << "   "
             << pos.hour << ":" << (pos.minute < 10 ? "0" : "") << pos.minute << "   "
             << pos.X << "   " << pos.Y << "   " << pos.Z << endl;
    }
    
    // 比较广播星历与精密星历
    cout << "\nComparing broadcast and precise ephemeris..." << endl;
    vector<PositionComparison> comparisons = CompareEphemeris(allPositions, precisePoints);
    
    if (!comparisons.empty()) {
        cout << "Found " << comparisons.size() << " matching points for comparison." << endl;
        
        // 保存比较结果
        SaveComparisonResults(comparisons, "ephemeris_comparison_results.txt");
        
        // 输出统计信息
        double sumDist = 0;
        double maxDist = 0;
        double minDist = 1e9;
        
        for (const auto& comp : comparisons) {
            sumDist += comp.distance;
            if (comp.distance > maxDist) maxDist = comp.distance;
            if (comp.distance < minDist) minDist = comp.distance;
        }
        
        double avgDist = sumDist / comparisons.size();
        
        cout << "\nComparison Statistics:" << endl;
        cout << "Average position difference: " << avgDist << " m" << endl;
        cout << "Maximum position difference: " << maxDist << " m" << endl;
        cout << "Minimum position difference: " << minDist << " m" << endl;
    } else {
        cout << "Warning: No matching points found for comparison!" << endl;
        cout << "Make sure the broadcast and precise ephemeris cover the same time period." << endl;
    }
    
    // 生成MATLAB绘图代码
    GenerateMatlabPlotCode(allPositions, "satellite_trajectories_3d.m");
    cout << "\nMATLAB plotting code generated: satellite_trajectories_3d.m" << endl;
    cout << "Run this file in MATLAB to generate the 3D trajectory plot." << endl;
    
    // 清理内存
    delete[] ephemerisBlocks;
    
    return 0;
}