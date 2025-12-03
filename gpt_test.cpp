#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

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
    double Toe, Cic, OMEGA, Cis;
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
        eph.OMEGA  = parseLineDouble(line4, 42, 19);
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

// ------------------- 输出到 TXT -------------------
void saveTxt(const string& filename, const vector<EphBlock>& data)
{
    ofstream fout(filename);

    for (const auto &e : data)
    {
        fout << fixed << setprecision(10);
        fout << e.sys << setw(2) << e.prn << " "
             << e.year << "-" << e.month << "-" << e.day << " "
             << e.hour << ":" << e.minute << ":" << e.second << "\n";

        fout << "a0 = " << e.a0 << "  a1 = " << e.a1 << "  a2 = " << e.a2 << "\n";
        fout << "e = " << e.e << "  sqrtA = " << e.sqrtA << "  Toe = " << e.Toe << "\n";
        fout << "i0 = " << e.i0 << "  omega = " << e.omega << "\n";
        fout << "---------------------------------------\n";
    }

    fout.close();
}

// ------------------- 主函数 -------------------
int main()
{
    string filename = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/brdm3350.19p";  // 你可以替换成任意 RINEX 文件

    vector<EphBlock> gps, bds;

    readNavFile(filename, gps, bds);

    saveTxt("gps_ephemeris.txt", gps);
    saveTxt("bds_ephemeris.txt", bds);

    cout << "读取完成！" << endl;
    cout << "GPS 星历数目: " << gps.size() << endl;
    cout << "BDS 星历数目: " << bds.size() << endl;

    return 0;
}
