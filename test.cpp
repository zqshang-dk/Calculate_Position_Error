#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// RINEX广播星历数据结构
struct BroadcastEphemeris {
    // 公共字段
    string satType;      // 卫星类型：G=GPS, C=BDS
    int prn;                  // 卫星PRN号
    int year, month, day, hour, minute;
    double second;
    
    // GPS/BDS公共的时钟参数
    double a0;               // 时钟偏差(s)
    double a1;               // 时钟漂移(s/s)
    double a2;               // 时钟漂移率(s/s?)
    
    // 轨道参数（适用于GPS和BDS）
    double Crs;              // 轨道半径正弦调和项振幅(m)
    double Delta_n;          // 平均运动差(rad/s)
    double M0;               // 平近点角(rad)
    
    double Cuc;              // 纬度幅角余弦调和项振幅(rad)
    double e;                // 轨道偏心率
    double Cus;              // 纬度幅角正弦调和项振幅(rad)
    double sqrtA;            // 轨道长半轴平方根(√m)
    
    double Toe;              // 星历参考时间(s of week)
    double Cic;              // 轨道倾角余弦调和项振幅(rad)
    double Omega0;           // 升交点赤经(rad)
    double Cis;              // 轨道倾角正弦调和项振幅(rad)
    
    double i0;               // 轨道倾角(rad)
    double Crc;              // 轨道半径余弦调和项振幅(m)
    double omega;            // 近地点幅角(rad)
    double OmegaDot;         // 升交点赤经变化率(rad/s)
    
    double IDOT;             // 轨道倾角变化率(rad/s)
    double L2Codes;          // L2上的码
    double GPSWeek;          // GPS周
    double L2PDataFlag;      // L2 P码数据标志
    
    double svAccuracy;       // 卫星精度(m)
    double svHealth;         // 卫星健康状况
    double TGD;              // 群延迟(s) - GPS使用，BDS可能为TGD1
    
    // BDS特有参数
    double TGD2;             // BDS群延迟2(s)
    double AODC;             // 时钟数据龄期
    double AODE;             // 星历数据龄期
    
    // 其他
    double IODC;             // 时钟数据IOD
    double IODE;             // 星历数据IOD
    double TransmissionTime; // 信息传输时间(s of week)
};

// 辅助函数：安全地转换字符串为double
double safeStod(const string& str) {
    if (str.empty() || str.find_first_not_of(' ') == string::npos) {
        return 0.0;
    }
    
    try {
        // 替换D为E（某些RINEX文件使用D表示指数）
        string processed = str;
        replace(processed.begin(), processed.end(), 'D', 'E');
        replace(processed.begin(), processed.end(), 'd', 'e');
        
        // 去除首尾空格
        size_t start = processed.find_first_not_of(" \t");
        size_t end = processed.find_last_not_of(" \t");
        
        if (start == string::npos || end == string::npos) {
            return 0.0;
        }
        
        return stod(processed.substr(start, end - start + 1));
    } catch (const exception&) {
        return 0.0;
    }
}

// 辅助函数：安全地转换字符串为int
int safeStoi(const string& str) {
    if (str.empty() || str.find_first_not_of(' ') == string::npos) {
        return 0;
    }
    
    try {
        // 去除首尾空格
        size_t start = str.find_first_not_of(" \t");
        size_t end = str.find_last_not_of(" \t");
        
        if (start == string::npos || end == string::npos) {
            return 0;
        }
        
        return stoi(str.substr(start, end - start + 1));
    } catch (const exception&) {
        return 0;
    }
}


// 读取广播星历文件
vector<BroadcastEphemeris> readRinexNav(const string& filename) {
    vector<BroadcastEphemeris> ephemerisList;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return ephemerisList;
    }
    
    string line;
    bool headerEnd = false;
    
    // 跳过文件头
    while (getline(file, line)) {
        // 检查是否到达头文件结尾
        if (line.find("END OF HEADER") != string::npos) {
            headerEnd = true;
            break;
        }
    }
    
    if (!headerEnd) {
        cerr << "Error: Could not find END OF HEADER" << endl;
        return ephemerisList;
    }
    
    // 读取星历数据
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        BroadcastEphemeris eph;
        
        try {
            // 第一行：卫星标识和时间信息
            if (line.length() < 23) continue;
            
            // 解析卫星类型和PRN号
            string satID = line.substr(0, 3);
            eph.satType = satID.substr(0, 1);  // G或C
            eph.prn = safeStoi(satID.substr(1, 2));
            
            // 解析时间
            eph.year = safeStoi(line.substr(4, 4));
            eph.month = safeStoi(line.substr(9, 2));
            eph.day = safeStoi(line.substr(12, 2));
            eph.hour = safeStoi(line.substr(15, 2));
            eph.minute = safeStoi(line.substr(18, 2));
            // 注意：秒可能有多位小数，我们取到下一个空格之前
            string secondStr = line.substr(21);
            size_t secEnd = secondStr.find(' ');
            if (secEnd != string::npos) {
                secondStr = secondStr.substr(0, secEnd);
            }
            eph.second = safeStod(line.substr(21, 2));
            
            // 时钟参数
            size_t startPos = 23;
            eph.a0 = stod(line.substr(startPos, 19));
            startPos += 19;
            eph.a1 = stod(line.substr(startPos, 19));
            startPos += 19;
            eph.a2 = stod(line.substr(startPos, 19));
            
            // 读取第二行
            if (!getline(file, line)) break;
            eph.Crs = stod(line.substr(23, 19));
            eph.Delta_n = stod(line.substr(42, 19));
            eph.M0 = stod(line.substr(61, 19));
            
            // 读取第三行
            if (!getline(file, line)) break;
            eph.Cuc = stod(line.substr(4, 19));
            eph.e = stod(line.substr(23, 19));
            eph.Cus = stod(line.substr(42, 19));
            eph.sqrtA = stod(line.substr(61, 19));
            
            // 读取第四行
            if (!getline(file, line)) break;
            eph.Toe = stod(line.substr(4, 19));
            eph.Cic = stod(line.substr(23, 19));
            eph.Omega0 = stod(line.substr(42, 19));
            eph.Cis = stod(line.substr(61, 19));
            
            // 读取第五行
            if (!getline(file, line)) break;
            eph.i0 = stod(line.substr(4, 19));
            eph.Crc = stod(line.substr(23, 19));
            eph.omega = stod(line.substr(42, 19));
            eph.OmegaDot = stod(line.substr(61, 19));
            
            // 读取第六行
            if (!getline(file, line)) break;
            eph.IDOT = stod(line.substr(4, 19));
            eph.L2Codes = stod(line.substr(23, 19));
            eph.GPSWeek = stod(line.substr(42, 19));
            eph.L2PDataFlag = stod(line.substr(61, 19));
            
            // 读取第七行
            if (!getline(file, line)) break;
            eph.svAccuracy = stod(line.substr(4, 19));
            eph.svHealth = stod(line.substr(23, 19));
            eph.TGD = stod(line.substr(42, 19));
            eph.IODC = stod(line.substr(61, 19));
            
            // 读取第八行
            if (!getline(file, line)) break;
            eph.TransmissionTime = stod(line.substr(4, 19));
            
            // 如果是BDS卫星，读取额外的参数
            if (eph.satType == "C") {
                eph.TGD2 = stod(line.substr(23, 19));
                eph.AODC = stod(line.substr(42, 19));
                eph.AODE = stod(line.substr(61, 19));
                eph.IODE = eph.AODE;  // BDS中AODE通常等于IODE
            } else {
                // GPS使用IODE字段
                if (line.length() > 61) {
                    eph.IODE = stod(line.substr(61, 19));
                }
            }
            
            ephemerisList.push_back(eph);
            
            // 跳过可能存在的额外空行
            while (file.peek() == '\n') {
                file.ignore(1);
            }
            
        } catch (const exception& e) {
            cerr << "Error parsing line: " << e.what() << endl;
            continue;
        }
    }
    
    file.close();
    return ephemerisList;
}

// 分离GPS和BDS星历数据
void separateGNSSData(const vector<BroadcastEphemeris>& allEph,
                      vector<BroadcastEphemeris>& gpsEph,
                      vector<BroadcastEphemeris>& bdsEph) {
    for (const auto& eph : allEph) {
        if (eph.satType == "G") {
            gpsEph.push_back(eph);
        } else if (eph.satType == "C") {
            bdsEph.push_back(eph);
        }
    }
}

// // 打印星历信息
// void printEphemerisInfo(const BroadcastEphemeris& eph) {
//     cout << fixed << setprecision(12);
//     cout << "\n=== " << eph.satType << setw(2) << setfill('0') 
//               << eph.prn << " Ephemeris ===" << endl;
//     cout << "Time: " << eph.year << "-" << eph.month << "-" << eph.day 
//               << " " << eph.hour << ":" << eph.minute << ":" << eph.second << endl;
//     cout << "Clock: a0=" << eph.a0 << ", a1=" << eph.a1 << ", a2=" << eph.a2 << endl;
//     cout << "Orbit: sqrtA=" << eph.sqrtA << ", e=" << eph.e << ", M0=" << eph.M0 << endl;
//     cout << "       Omega0=" << eph.Omega0 << ", i0=" << eph.i0 << ", omega=" << eph.omega << endl;
//     cout << "Toe=" << eph.Toe << ", GPSWeek=" << eph.GPSWeek << endl;
    
//     if (eph.satType == "C") {
//         cout << "BDS Specific: TGD1=" << eph.TGD << ", TGD2=" << eph.TGD2 
//                   << ", AODC=" << eph.AODC << ", AODE=" << eph.AODE << endl;
//     } else {
//         cout << "GPS Specific: TGD=" << eph.TGD << ", IODC=" << eph.IODC 
//                   << ", IODE=" << eph.IODE << endl;
//     }
// }

// 将GPS星历信息写入文件
void writeGPSEphemerisToFile(const vector<BroadcastEphemeris>& gpsEph, const string& filename) {
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error: Could not open GPS output file " << filename << endl;
        return;
    }
    
    outFile << fixed << setprecision(12);
    outFile << "=== GPS Ephemeris Data ===\n";
    outFile << "Total records: " << gpsEph.size() << "\n\n";
    
    for (const auto& eph : gpsEph) {
        outFile << "=== G" << setw(2) << setfill('0') << eph.prn << " Ephemeris ===\n";
        outFile << "Time: " << eph.year << "-" << eph.month << "-" << eph.day 
                << " " << eph.hour << ":" << eph.minute << ":" << eph.second << "\n";
        outFile << "Clock: a0=" << eph.a0 << ", a1=" << eph.a1 << ", a2=" << eph.a2 << "\n";
        outFile << "Orbit: sqrtA=" << eph.sqrtA << ", e=" << eph.e << ", M0=" << eph.M0 << "\n";
        outFile << "       Omega0=" << eph.Omega0 << ", i0=" << eph.i0 << ", omega=" << eph.omega << "\n";
        outFile << "       Delta_n=" << eph.Delta_n << ", OmegaDot=" << eph.OmegaDot << ", IDOT=" << eph.IDOT << "\n";
        outFile << "       Crs=" << eph.Crs << ", Crc=" << eph.Crc << ", Cuc=" << eph.Cuc << ", Cus=" << eph.Cus << "\n";
        outFile << "       Cic=" << eph.Cic << ", Cis=" << eph.Cis << "\n";
        outFile << "Toe=" << eph.Toe << ", GPSWeek=" << eph.GPSWeek << "\n";
        outFile << "GPS Specific: TGD=" << eph.TGD << ", IODC=" << eph.IODC 
                << ", IODE=" << eph.IODE << "\n";
        outFile << "svAccuracy=" << eph.svAccuracy << ", svHealth=" << eph.svHealth << "\n";
        outFile << "TransmissionTime=" << eph.TransmissionTime << "\n";
        outFile << "=====================================\n\n";
    }
    
    outFile.close();
    cout << "GPS ephemeris data written to " << filename << endl;
}

// 将BDS星历信息写入文件
void writeBDSEphemerisToFile(const vector<BroadcastEphemeris>& bdsEph, const string& filename) {
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error: Could not open BDS output file " << filename << endl;
        return;
    }
    
    outFile << fixed << setprecision(12);
    outFile << "=== BDS (BeiDou) Ephemeris Data ===\n";
    outFile << "Total records: " << bdsEph.size() << "\n\n";
    
    for (const auto& eph : bdsEph) {
        outFile << "=== C" << setw(2) << setfill('0') << eph.prn << " Ephemeris ===\n";
        outFile << "Time: " << eph.year << "-" << eph.month << "-" << eph.day 
                << " " << eph.hour << ":" << eph.minute << ":" << eph.second << "\n";
        outFile << "Clock: a0=" << eph.a0 << ", a1=" << eph.a1 << ", a2=" << eph.a2 << "\n";
        outFile << "Orbit: sqrtA=" << eph.sqrtA << ", e=" << eph.e << ", M0=" << eph.M0 << "\n";
        outFile << "       Omega0=" << eph.Omega0 << ", i0=" << eph.i0 << ", omega=" << eph.omega << "\n";
        outFile << "       Delta_n=" << eph.Delta_n << ", OmegaDot=" << eph.OmegaDot << ", IDOT=" << eph.IDOT << "\n";
        outFile << "       Crs=" << eph.Crs << ", Crc=" << eph.Crc << ", Cuc=" << eph.Cuc << ", Cus=" << eph.Cus << "\n";
        outFile << "       Cic=" << eph.Cic << ", Cis=" << eph.Cis << "\n";
        outFile << "Toe=" << eph.Toe << ", GPSWeek=" << eph.GPSWeek << "\n";
        outFile << "BDS Specific: TGD1=" << eph.TGD << ", TGD2=" << eph.TGD2 
                << ", AODC=" << eph.AODC << ", AODE=" << eph.AODE << "\n";
        outFile << "svAccuracy=" << eph.svAccuracy << ", svHealth=" << eph.svHealth << "\n";
        outFile << "TransmissionTime=" << eph.TransmissionTime << "\n";
        outFile << "=====================================\n\n";
    }
    
    outFile.close();
    cout << "BDS ephemeris data written to " << filename << endl;
}

// 示例使用
int main() {
    string filename = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/brdm3350.19p";
    
    // 读取所有星历数据
    vector<BroadcastEphemeris> allEphemeris = readRinexNav(filename);
    
    if (allEphemeris.empty()) {
        cout << "No ephemeris data found." << endl;
        return 1;
    }
    
    cout << "Total ephemeris records: " << allEphemeris.size() << endl;
    
    // 分离GPS和BDS数据
    vector<BroadcastEphemeris> gpsEphemeris;
    vector<BroadcastEphemeris> bdsEphemeris;
    
    separateGNSSData(allEphemeris, gpsEphemeris, bdsEphemeris);
    
    cout << "\nGPS satellites: " << gpsEphemeris.size() << endl;
    cout << "BDS satellites: " << bdsEphemeris.size() << endl;
    
    writeGPSEphemerisToFile(gpsEphemeris, "gps_ephemeris.txt");
    writeBDSEphemerisToFile(bdsEphemeris, "bds_ephemeris.txt");
    
    
    // 保存分离的数据以便后续处理
    cout << "\nData separation complete. You can now use gpsEphemeris and bdsEphemeris vectors for position calculation." << endl;
    
    return 0;
}