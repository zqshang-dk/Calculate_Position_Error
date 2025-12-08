#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense> // 依赖于 Eigen 库
#include <set>
#include <map>
#include <algorithm>

// 定义一个结构体来存储卫星位置数据
struct SatellitePosition {
    int year, month, day, hour, minute;
    double second;
    std::string prn;
    Eigen::Vector3d position; // X, Y, Z
    double clock;           // 钟差 (m)
};

// 重新定义 Key 的类型：{从2019/12/1 00:00开始的总分钟数, PRN}
using PositionKey = std::pair<int, std::string>;

// 用于存储精密星历数据的映射
using PreciseEphemerisMap = std::map<PositionKey, SatellitePosition>;
// 用于存储广播星历数据的映射
using BroadcastEphemerisMap = std::map<PositionKey, SatellitePosition>;

/**
 * @brief 计算时间点相对于 2019-12-01 00:00:00 的总分钟数
 * @param pos 卫星位置数据结构体
 * @return 总分钟数
 */
int CalculateTotalMinutes(const SatellitePosition& pos) {
    // 2019 12 1 0 0 0.00000000 作为起始点
    // 简单地计算从 00:00 开始的分钟数，因为我们只处理 1 天的数据
    // 如果处理多天，需要加上 (pos.day - 1) * 24 * 60
    return pos.hour * 60 + pos.minute;
}

/**
 * @brief 读取 SP3 精密星历文件
 * @param filename SP3 文件路径
 * @return 存储精密星历数据的映射
 */
// PreciseEphemerisMap ReadSP3(const std::string& filename) {
//     PreciseEphemerisMap precise_map;
//     std::ifstream file(filename);
//     std::string line;
    
//     if (!file.is_open()) {
//         std::cerr << "Error: Unable to open SP3 file: " << filename << std::endl;
//         return precise_map;
//     }

//     const double KM_TO_M = 1000.0;

//     while (std::getline(file, line)) {
//         if (line.length() < 3) continue;

//         // 读取时间标记 (Epoch Flag)
//         if (line[0] == '*') {
//             int year, month, day, hour, minute;
//             double second;
            
//             // 假设 SP3 文件时间格式为：* YYYY MM DD HH MM SS.sssss...
//             if (line.length() >= 25 && isdigit(line[1])) {
//                 std::stringstream ss(line.substr(2));
//                 // 确保时间读取的正确性
//                 if (!(ss >> year >> month >> day >> hour >> minute >> second)) {
//                     continue; // 读取失败则跳过
//                 }

//                 // 只处理 15min 间隔的数据，并忽略秒数非 0 的情况
//                 if (minute % 15 != 0 || std::abs(second) > 1e-6) {
//                     continue;
//                 }

//                 SatellitePosition time_key = {year, month, day, hour, minute, second};
//                 int current_total_minutes = CalculateTotalMinutes(time_key);

//                 // 在下一个循环中读取该时间戳下的卫星位置
//                 while (std::getline(file, line) && line[0] == 'P') {
//                     if (line.length() < 50) continue;

//                     std::string prn = line.substr(1, 3);
//                     double x, y, z, clk;
//                     std::stringstream ss_pos(line.substr(4));
                    
//                     // 只处理 G 和 C 卫星 (GPS 和 BDS)
//                     if (prn[0] == 'G' || prn[0] == 'C') {
//                         if (!(ss_pos >> x >> y >> z >> clk)) {
//                             continue; // 位置数据读取失败则跳过
//                         }

//                         SatellitePosition pos = time_key;
//                         pos.prn = prn;
//                         pos.position = Eigen::Vector3d(x * KM_TO_M, y * KM_TO_M, z * KM_TO_M);
//                         pos.clock = clk * KM_TO_M * 1e-6; // 转换为 m

//                         // 修正后的 insert 语句
//                         PositionKey key = {current_total_minutes, pos.prn};
//                         precise_map.insert({key, pos}); // 使用 std::pair<Key, Value> 结构
//                     }
//                 }
//             }
//         }
//     }
//     file.close();
//     std::cout << "Successfully read " << precise_map.size() << " position records from SP3 file." << std::endl;
//     return precise_map;
// }
/**
 * @brief 读取 SP3 精密星历文件 (改进版本：更稳健地跳过头部)
 * @param filename SP3 文件路径
 * @return 存储精密星历数据的映射
 */
PreciseEphemerisMap ReadSP3(const std::string& filename) {
    PreciseEphemerisMap precise_map;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open SP3 file: " << filename << std::endl;
        return precise_map;
    }

    const double KM_TO_M = 1000.0;
    bool in_header = true; // 增加一个标志，用来跳过所有头部信息

    while (std::getline(file, line)) {
        if (line.length() < 3) continue;

        // 1. 跳过头部注释行
        if (in_header) {
            char first_char = line[0];
            if (first_char == '#' || first_char == '+' || first_char == '%' || first_char == '/' || first_char == ' ') {
                continue; // 忽略头部行
            }
            // 找到第一个 * 行，表示头部结束，开始读取数据
            if (first_char == '*') {
                in_header = false;
            } else {
                continue; // 依然在头部，但格式不符，跳过
            }
        }
        
        // 2. 读取时间标记 (Epoch Flag)
        if (line[0] == '*') {
            int year, month, day, hour, minute;
            double second;
            
            // SP3 标准时间格式：* YYYY MM DD HH MM SS.sssss...
            if (line.length() >= 25) {
                // 确保是从第 3 个字符开始读取，跳过 * 和空格
                std::stringstream ss(line.substr(2));
                
                if (!(ss >> year >> month >> day >> hour >> minute >> second)) {
                    continue; // 时间读取失败则跳过
                }

                // 只处理 15min 间隔的数据，且秒数为 0
                if (minute % 15 != 0 || std::abs(second) > 1e-6) {
                    continue;
                }

                SatellitePosition time_key = {year, month, day, hour, minute, second};
                int current_total_minutes = CalculateTotalMinutes(time_key);

                // 在下一个循环中读取该时间戳下的卫星位置
                while (std::getline(file, line) && line.length() > 0 && line[0] == 'P') {
                    if (line.length() < 50) continue;

                    std::string prn = line.substr(1, 3);
                    double x, y, z, clk;
                    // 位置数据从第 5 个字符开始读取 (P G01 )
                    std::stringstream ss_pos(line.substr(4));
                    
                    if (!(ss_pos >> x >> y >> z >> clk)) {
                        continue; // 位置数据读取失败则跳过
                    }
                    
                    // 只处理 G 和 C 卫星 (GPS 和 BDS)
                    if (prn[0] == 'G' || prn[0] == 'C') {
                        SatellitePosition pos = time_key;
                        pos.prn = prn;
                        pos.position = Eigen::Vector3d(x * KM_TO_M, y * KM_TO_M, z * KM_TO_M);
                        pos.clock = clk * KM_TO_M * 1e-6; // 转换为 m

                        PositionKey key = {current_total_minutes, pos.prn};
                        precise_map.insert({key, pos}); 
                    }
                }
            }
        }
    }
    file.close();
    std::cout << "Successfully read " << precise_map.size() << " position records from SP3 file." << std::endl;
    return precise_map;
}

/**
 * @brief 读取广播星历计算的位置文件
 * @param filename 文件路径
 * @return 存储广播星历计算位置数据的映射
 */
BroadcastEphemerisMap ReadBroadcastEphemeris(const std::string& filename) {
    BroadcastEphemerisMap broadcast_map;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open broadcast ephemeris file: " << filename << std::endl;
        return broadcast_map;
    }

    SatellitePosition current_time_key;
    
    while (std::getline(file, line)) {
        if (line.length() == 0) continue;

        if (line[0] == '*') {
            // 时间行: * 2019 12 1 0 0 0.00000000 900.00000000
            std::stringstream ss(line.substr(1));
            // 确保时间读取的正确性
            if (!(ss >> current_time_key.year >> current_time_key.month >> current_time_key.day 
               >> current_time_key.hour >> current_time_key.minute >> current_time_key.second)) {
                   continue;
            }
            
            current_time_key.second = 0.0;
            
            continue;
        }

        // 卫星位置行: G02 -15007423.16652131 14399954.52302857 17077873.16419729 4397.88031535
        if (line[0] == 'G' || line[0] == 'C') {
            std::stringstream ss(line);
            SatellitePosition pos = current_time_key;
            double x, y, z, clk;

            if (!(ss >> pos.prn >> x >> y >> z >> clk)) {
                continue; // 数据读取失败则跳过
            }
            
            pos.position = Eigen::Vector3d(x, y, z); // 单位是 m
            pos.clock = clk;                         // 单位是 m
            
            // 修正后的 insert 语句
            int current_total_minutes = CalculateTotalMinutes(pos);
            PositionKey key = {current_total_minutes, pos.prn};
            broadcast_map.insert({key, pos});
        }
    }
    file.close();
    std::cout << "Successfully read " << broadcast_map.size() << " position records from broadcast ephemeris file." << std::endl;
    return broadcast_map;
}


/**
 * @brief 对比精密星历和广播星历的位置，计算误差并输出到文件
 */
void CalculateError(const PreciseEphemerisMap& precise_map, const BroadcastEphemerisMap& broadcast_map, char system, const std::string& output_filename) {
    std::ofstream outfile(output_filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open output file: " << output_filename << std::endl;
        return;
    }

    // 写入文件头
    outfile << std::fixed << std::setprecision(4);
    outfile << "# Time (min)   PRN   dX (m)   dY (m)   dZ (m)   3D Error (m)   dClk (m)" << std::endl;
    outfile << "# ----------------------------------------------------------------------" << std::endl;
    
    std::set<std::string> target_prns;
    int error_count = 0;

    // 遍历广播星历数据，以广播星历的时间/PRN为基准进行查找
    for (const auto& pair : broadcast_map) {
        const PositionKey& key = pair.first;
        const SatellitePosition& broadcast_pos = pair.second;

        if (broadcast_pos.prn[0] != system) continue;

        // 查找对应时间的精密星历数据
        if (precise_map.count(key)) {
            const SatellitePosition& precise_pos = precise_map.at(key);

            // 1. 计算位置误差 (广播星历 - 精密星历)
            Eigen::Vector3d error_vector = broadcast_pos.position - precise_pos.position;

            // 2. 计算三维误差 (3D Error)
            double error_3d = error_vector.norm();

            // 3. 计算钟差误差 (dClk)
            double error_clk = broadcast_pos.clock - precise_pos.clock;

            // 输出结果
            outfile << std::setw(12) << key.first // 总分钟数
                    << std::setw(6) << key.second // PRN
                    << std::setw(10) << error_vector(0)
                    << std::setw(10) << error_vector(1)
                    << std::setw(10) << error_vector(2)
                    << std::setw(15) << error_3d
                    << std::setw(12) << error_clk
                    << std::endl;
            
            target_prns.insert(key.second);
            error_count++;
        }
    }

    outfile.close();
    std::cout << "Successfully calculated and saved " << system << " system errors to " << output_filename << ". Total matched records: " << error_count << std::endl;
}

/**
 * @brief 主逻辑函数，执行数据读取和误差计算
 */
void CalculatePositionError() {
    // 假设您已将文件放在程序运行目录下
    std::string sp3_file = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/WUM0MGXFIN_20193350000_01D_15M_ORB.SP3";
    std::string gps_broadcast_file = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/gps_position.txt";
    std::string bds_broadcast_file = "E:/STUDY/Sophomore1/卫星导航原理/Calculate_Position_Error/bds_position.txt";
    std::string gps_error_output = "gps_orbit_error.txt";
    std::string bds_error_output = "bds_orbit_error.txt";

    // --- 1. 读取精密星历文件（SP3） ---
    PreciseEphemerisMap precise_ephemeris = ReadSP3(sp3_file);

    // --- 2. 读取广播星历计算的位置文件 ---
    BroadcastEphemerisMap gps_broadcast_ephemeris = ReadBroadcastEphemeris(gps_broadcast_file);
    BroadcastEphemerisMap bds_broadcast_ephemeris = ReadBroadcastEphemeris(bds_broadcast_file);

    // --- 3. 计算 GPS 误差 ---
    CalculateError(precise_ephemeris, gps_broadcast_ephemeris, 'G', gps_error_output);

    // --- 4. 计算 BDS 误差 ---
    CalculateError(precise_ephemeris, bds_broadcast_ephemeris, 'C', bds_error_output);
}

// 主函数
int main() {
    // 设置输出流的浮点数格式
    std::cout << std::fixed << std::setprecision(2);
    
    CalculatePositionError();

    return 0;
}