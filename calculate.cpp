#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

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
    double time;      // 时间 (秒)
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

// GPS卫星位置计算
SatellitePosition CalculateGpsSatellitePosition(const EPHEMERISBLOCK& eph, double transmitTime) {
    SatellitePosition pos;
    pos.PRN = eph.PRN;
    pos.satSystem = 'G';
    pos.time = transmitTime;
    
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
SatellitePosition CalculateBdsSatellitePosition(const EPHEMERISBLOCK& eph, double transmitTime) {
    SatellitePosition pos;
    pos.PRN = eph.PRN;
    pos.satSystem = 'C';
    pos.time = transmitTime;
    
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

// 生成MATLAB绘图代码
void GenerateMatlabPlotCode(const vector<SatellitePosition>& positions, const char* filename) {
    ofstream matlabFile(filename);
    
    matlabFile << "% MATLAB code for satellite trajectory plotting" << endl;
    matlabFile << "clear; close all; clc;" << endl << endl;
    
    // 分离GPS和BDS卫星数据
    vector<int> gpsPRNs, bdsPRNs;
    for (const auto& pos : positions) {
        if (pos.satSystem == 'G') {
            if (find(gpsPRNs.begin(), gpsPRNs.end(), pos.PRN) == gpsPRNs.end()) {
                gpsPRNs.push_back(pos.PRN);
            }
        } else if (pos.satSystem == 'C') {
            if (find(bdsPRNs.begin(), bdsPRNs.end(), pos.PRN) == bdsPRNs.end()) {
                bdsPRNs.push_back(pos.PRN);
            }
        }
    }
    
    // 组织数据
    matlabFile << "% Satellite positions data" << endl;
    matlabFile << "times = [";
    for (const auto& pos : positions) {
        matlabFile << pos.time << " ";
    }
    matlabFile << "];" << endl << endl;
    
    // GPS卫星数据
    for (int prn : gpsPRNs) {
        matlabFile << "gps" << prn << "_X = [";
        for (const auto& pos : positions) {
            if (pos.satSystem == 'G' && pos.PRN == prn) {
                matlabFile << pos.X << " ";
            }
        }
        matlabFile << "];" << endl;
        
        matlabFile << "gps" << prn << "_Y = [";
        for (const auto& pos : positions) {
            if (pos.satSystem == 'G' && pos.PRN == prn) {
                matlabFile << pos.Y << " ";
            }
        }
        matlabFile << "];" << endl;
        
        matlabFile << "gps" << prn << "_Z = [";
        for (const auto& pos : positions) {
            if (pos.satSystem == 'G' && pos.PRN == prn) {
                matlabFile << pos.Z << " ";
            }
        }
        matlabFile << "];" << endl << endl;
    }
    
    // BDS卫星数据
    for (int prn : bdsPRNs) {
        matlabFile << "bds" << prn << "_X = [";
        for (const auto& pos : positions) {
            if (pos.satSystem == 'C' && pos.PRN == prn) {
                matlabFile << pos.X << " ";
            }
        }
        matlabFile << "];" << endl;
        
        matlabFile << "bds" << prn << "_Y = [";
        for (const auto& pos : positions) {
            if (pos.satSystem == 'C' && pos.PRN == prn) {
                matlabFile << pos.Y << " ";
            }
        }
        matlabFile << "];" << endl;
        
        matlabFile << "bds" << prn << "_Z = [";
        for (const auto& pos : positions) {
            if (pos.satSystem == 'C' && pos.PRN == prn) {
                matlabFile << pos.Z << " ";
            }
        }
        matlabFile << "];" << endl << endl;
    }
    
    // 绘图代码
    matlabFile << "% Plot satellite trajectories" << endl;
    matlabFile << "figure('Position', [100, 100, 1200, 800]);" << endl;
    
    // 3D轨迹图
    matlabFile << "subplot(2,2,1);" << endl;
    matlabFile << "hold on; grid on;" << endl;
    
    for (int prn : gpsPRNs) {
        matlabFile << "plot3(gps" << prn << "_X, gps" << prn << "_Y, gps" << prn << "_Z, 'r-', 'LineWidth', 1.5);" << endl;
    }
    
    for (int prn : bdsPRNs) {
        matlabFile << "plot3(bds" << prn << "_X, bds" << prn << "_Y, bds" << prn << "_Z, 'b-', 'LineWidth', 1.5);" << endl;
    }
    
    matlabFile << "xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');" << endl;
    matlabFile << "title('Satellite Trajectories in ECEF');" << endl;
    matlabFile << "legend('GPS', 'BDS');" << endl;
    matlabFile << "axis equal;" << endl;
    
    // XY平面投影
    matlabFile << "subplot(2,2,2);" << endl;
    matlabFile << "hold on; grid on;" << endl;
    
    for (int prn : gpsPRNs) {
        matlabFile << "plot(gps" << prn << "_X, gps" << prn << "_Y, 'r-', 'LineWidth', 1.5);" << endl;
    }
    
    for (int prn : bdsPRNs) {
        matlabFile << "plot(bds" << prn << "_X, bds" << prn << "_Y, 'b-', 'LineWidth', 1.5);" << endl;
    }
    
    matlabFile << "xlabel('X (m)'); ylabel('Y (m)');" << endl;
    matlabFile << "title('XY Plane Projection');" << endl;
    matlabFile << "axis equal;" << endl;
    
    matlabFile << "saveas(gcf, 'satellite_trajectories.png');" << endl;
    matlabFile << "disp('Plot saved as satellite_trajectories.png');" << endl;
    
    matlabFile.close();
}

int main() {
    // 读取广播星历文件
    const char* broadcastFile = "brdm3350.19p";
    int ephemerisBlockNum = 0;
    EPHEMERISBLOCK* ephemerisBlocks = nullptr;
    
    cout << "Reading broadcast ephemeris file..." << endl;
    if (ReadBrodcastEphemeris(broadcastFile, ephemerisBlockNum, ephemerisBlocks) != 0) {
        cerr << "Failed to read broadcast ephemeris file!" << endl;
        return -1;
    }
    cout << "Successfully read " << ephemerisBlockNum << " ephemeris blocks." << endl;
    
    // 计算当天任意时刻的卫星位置
    vector<SatellitePosition> allPositions;
    int calculationDate[] = {2019, 12, 1};  // 年月日
    
    // 计算从00:00:00到23:59:59，每隔5分钟计算一次
    for (int hour = 0; hour < 24; hour++) {
        for (int minute = 0; minute < 60; minute += 5) {
            double transmitTime = hour * 3600.0 + minute * 60.0;
            
            for (int i = 0; i < ephemerisBlockNum; i++) {
                SatellitePosition pos;
                
                if (ephemerisBlocks[i].satSystem == 'G') {
                    pos = CalculateGpsSatellitePosition(ephemerisBlocks[i], transmitTime);
                } else if (ephemerisBlocks[i].satSystem == 'C') {
                    pos = CalculateBdsSatellitePosition(ephemerisBlocks[i], transmitTime);
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
    cout << "PRN Sys Time(s)         X(m)           Y(m)           Z(m)" << endl;
    cout.precision(3);
    cout << fixed;
    
    for (int i = 0; i < min(10, (int)allPositions.size()); i++) {
        const auto& pos = allPositions[i];
        cout << pos.PRN << "   " << pos.satSystem << "   " 
             << pos.time << "   " << pos.X << "   " << pos.Y << "   " << pos.Z << endl;
    }
    
    // 生成MATLAB绘图代码
    GenerateMatlabPlotCode(allPositions, "satellite_trajectories.m");
    cout << "\nMATLAB plotting code generated: satellite_trajectories.m" << endl;
    
    // 清理内存
    delete[] ephemerisBlocks;
    
    return 0;
}