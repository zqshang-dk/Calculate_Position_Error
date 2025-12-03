%% 完整的广播星历vs精密星历对比代码
% 适配你的SP3格式（历元头：年月日时分秒）和广播星历格式
% 运行前仅需修改：文件路径部分

close all; clear; clc;

%% ==================== 1. 自定义函数定义 ====================
% --------------------- 读取广播星历文件（适配你的格式） ---------------------
function brdc = readBrdcSP3(brdc_file)
    % 输入：brdc_file - 广播星历文件名（如gps_position.txt）
    % 输出：brdc - 表，包含PRN、t(当日累计秒)、X/Y/Z(米)
    
    fid = fopen(brdc_file, 'r');
    if fid == -1
        error('无法打开广播星历文件：%s，请检查路径是否正确！', brdc_file);
    end
    
    % 初始化存储变量
    PRN = [];       % 卫星PRN号（数字，如1/2/3）
    t = [];         % 当日累计秒（0~86400）
    X = []; Y = []; Z = [];  % 卫星坐标（米）
    current_t = 0;  % 当前历元的累计秒
    
    % 逐行读取文件
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line) || strncmp(line, '#', 1)
            continue;  % 跳过空行/注释行
        end
        
        % 识别历元头行（格式：* 2019 12 1 0 0 0.00000000 900.00000000）
        if strncmp(line, '*', 1)
            parts = strsplit(strtrim(line));
            if length(parts) >= 7
                % 解析时分秒，计算当日累计秒
                hour = str2double(parts{5});
                minute = str2double(parts{6});
                sec = str2double(parts{7});
                current_t = hour*3600 + minute*60 + sec;
            end
            continue;
        end
        
        % 识别卫星数据行（如G02 -15007423.16652131 14399954.52302857...）
        parts = strsplit(strtrim(line));
        if length(parts) < 5
            continue;  % 跳过格式错误行
        end
        
        % 解析PRN（如G02→2，C05→5）
        prn_str = parts{1};
        prn_num = str2double(prn_str(2:end));
        if isnan(prn_num)
            continue;
        end
        
        % 解析坐标（米）
        x = str2double(parts{2});
        y = str2double(parts{3});
        z = str2double(parts{4});
        if any(isnan([x,y,z]))
            continue;
        end
        
        % 存储数据
        PRN = [PRN; prn_num];
        t = [t; current_t];
        X = [X; x];
        Y = [Y; y];
        Z = [Z; z];
    end
    fclose(fid);
    
    % 转换为Table格式（方便后续处理）
    if ~isempty(PRN)
        brdc = table(PRN, t, X, Y, Z, 'VariableNames', {'PRN','t','X','Y','Z'});
    else
        error('广播星历文件%s中无有效数据！', brdc_file);
    end
end

% --------------------- 读取精密星历SP3文件（修正版） ---------------------
function prec = readSP3(sp3_file)
    % 输入：sp3_file - 精密星历SP3文件名
    % 输出：prec - 表，包含PRN、t(当日累计秒)、X/Y/Z(米)
    
    fid = fopen(sp3_file, 'r');
    if fid == -1
        error('无法打开精密星历文件：%s，请检查路径是否正确！', sp3_file);
    end
    
    % 初始化存储变量
    PRN = [];       % 卫星PRN号（数字）
    t = [];         % 当日累计秒
    X = []; Y = []; Z = [];  % 卫星坐标（米）
    current_t = 0;  % 当前历元的累计秒
    line_count = 0;
    
    % 逐行读取文件
    while ~feof(fid)
        line = fgetl(fid);
        line_count = line_count + 1;
        if isempty(line)
            continue;
        end
        
        % 识别历元头行（格式：* 2019 12  1  0  0  0.00000000）
        if strncmp(line, '*', 1)
            parts = strsplit(strtrim(line));
            if length(parts) >= 7
                % 解析时分秒，计算当日累计秒
                hour = str2double(parts{4});
                minute = str2double(parts{5});
                sec = str2double(parts{6});
                current_t = hour*3600 + minute*60 + sec;
            end
            continue;
        end
        
        % 识别卫星数据行（格式以P开头，如：PG01  -8420.638553 -13775.264971 -21352.655083）
        if strncmp(line, 'P', 1) && length(line) >= 60
            try
                % 解析卫星标识（如PG01、PC05等）
                sat_id = strtrim(line(2:5));  % 提取PG01等
                if length(sat_id) < 2
                    continue;
                end
                
                % 解析系统类型和PRN号
                sys_char = sat_id(1);  % 系统标识：G=GPS, C=BDS, R=GLONASS, E=Galileo
                prn_str = sat_id(2:end);
                prn_num = str2double(prn_str);
                
                if isnan(prn_num)
                    continue;
                end
                
                % 解析坐标（SP3单位是km，转换为米）
                x_str = strtrim(line(5:20));
                y_str = strtrim(line(21:36));
                z_str = strtrim(line(37:52));
                
                x = str2double(x_str) * 1000;  % km转m
                y = str2double(y_str) * 1000;
                z = str2double(z_str) * 1000;
                
                if any(isnan([x,y,z]))
                    continue;
                end
                
                % 只保留GPS和BDS数据
                if strcmp(sys_char, 'G') || strcmp(sys_char, 'C')
                    PRN = [PRN; prn_num];
                    t = [t; current_t];
                    X = [X; x];
                    Y = [Y; y];
                    Z = [Z; z];
                end
            catch
                fprintf('警告：第%d行解析失败：%s\n', line_count, line);
                continue;
            end
        end
    end
    fclose(fid);
    
    % 转换为Table格式
    if ~isempty(PRN)
        prec = table(PRN, t, X, Y, Z, 'VariableNames', {'PRN','t','X','Y','Z'});
        fprintf('成功读取%d个精密星历数据点\n', length(PRN));
    else
        % 显示文件前几行帮助调试
        fid = fopen(sp3_file, 'r');
        fprintf('SP3文件前10行内容：\n');
        for i = 1:10
            if ~feof(fid)
                line = fgetl(fid);
                fprintf('行%d: %s\n', i, line);
            end
        end
        fclose(fid);
        error('精密星历文件%s中无有效GPS/BDS数据！请检查文件格式。', sp3_file);
    end
end

%% ==================== 2. 主逻辑：广播星历vs精密星历对比 ====================
% --------------------- 配置文件路径（仅需修改这里！） ---------------------
gps_brdc_file = 'gps_position.txt';    % GPS广播星历文件
bds_brdc_file = 'bds_position.txt';    % BDS广播星历文件
prec_sp3_file = 'WUM0MGXFIN_20193350000_01D_15M_ORB.SP3';  % 精密星历SP3文件

% --------------------- 读取数据 ---------------------
fprintf('正在读取广播星历数据...\n');
brdc_gps = readBrdcSP3(gps_brdc_file);
% brdc_bds = readBrdcSP3(bds_brdc_file);  % 先注释BDS，测试GPS正常后再放开

fprintf('正在读取精密星历数据...\n');
prec = readSP3(prec_sp3_file);

% --------------------- 初始化参数 ---------------------
systems = {'G'};       % 先测试GPS，后续可加'C'
sys_names = {'GPS'};   % 先测试GPS，后续可加'BDS'
all_rms = NaN(length(systems),1);        % 存储各系统RMS误差
time_tol = 15;             % 时间匹配容差（秒）

% --------------------- 创建绘图窗口 ---------------------
figure('Position',[100 100 1100 700], 'Color','w');
tiledlayout(length(systems),1,'TileSpacing','compact','Padding','compact');

% --------------------- 逐系统对比 ---------------------
for si = 1:length(systems)
    sys = systems{si};
    sys_name = sys_names{si};
    
    % 选择对应广播星历
    if strcmp(sys, 'G')
        brdc = brdc_gps;
    elseif strcmp(sys, 'C')
        brdc = brdc_bds;
    else
        error('不支持的系统类型：%s', sys);
    end
    
    % 提取唯一PRN和历元
    prns = unique(brdc.PRN);
    prns = sort(prns);
    epochs = unique(brdc.t);
    epochs = sort(epochs);
    
    % 初始化误差矩阵（行=历元，列=PRN）
    err_mat = NaN(length(epochs), length(prns));
    prn2col = containers.Map(prns, 1:length(prns));
    
    % 遍历每个历元计算误差
    fprintf('正在计算%s系统轨道误差...\n', sys_name);
    for e_idx = 1:length(epochs)
        t_epoch = epochs(e_idx);
        % 提取该历元的广播星历数据
        brdc_idx = brdc.t == t_epoch;
        brdc_epoch = brdc(brdc_idx, :);
        
        % 遍历该历元的所有卫星
        for i = 1:height(brdc_epoch)
            prn = brdc_epoch.PRN(i);
            if ~isKey(prn2col, prn)
                continue;
            end
            col = prn2col(prn);
            
            % 广播星历坐标
            Xb = [brdc_epoch.X(i), brdc_epoch.Y(i), brdc_epoch.Z(i)];
            
            % 匹配精密星历（同PRN + 时间差<容差）
            prec_idx = prec.PRN == prn & abs(prec.t - t_epoch) < time_tol;
            if any(prec_idx)
                % 精密星历坐标
                Xp = [prec.X(prec_idx), prec.Y(prec_idx), prec.Z(prec_idx)];
                % 计算3D轨道误差（米）
                err_mat(e_idx, col) = norm(Xb - Xp);
            end
        end
    end
    
    % --------------------- 绘制误差时序图 ---------------------
    nexttile;
    hold on; grid on; box on;
    colors = lines(length(prns));  % 自动生成区分度高的颜色
    
    % 逐卫星绘制
    for k = 1:length(prns)
        prn = prns(k);
        plot(epochs/3600, err_mat(:,k), ...
             'Color', colors(k,:), ...
             'LineWidth', 1.2, ...
             'DisplayName', sprintf('%s%02d', sys, prn));
    end
    
    % 图形美化
    xlabel('Epoch (Hour of Day)','FontSize',10);
    ylabel('3D Orbit Error (m)','FontSize',10);
    title(sprintf('%s Broadcast vs Precise Ephemeris Error', sys_name), 'FontSize',11);
    legend('show','Location','eastoutside','FontSize',8,'NumColumns',1);
    ylim([0, max(err_mat(:), [], 'omitNaN')*1.1]);  % 自适应Y轴
    set(gca, 'FontSize', 9);
    
    % --------------------- 计算RMS误差 ---------------------
    valid_err = err_mat(~isnan(err_mat));
    if ~isempty(valid_err)
        rms = sqrt(mean(valid_err.^2));
        all_rms(si) = rms;
        fprintf('%s 总体3D轨道RMS误差：%.3f 米\n', sys_name, rms);
    else
        fprintf('%s 无有效匹配的历元数据！\n', sys_name);
    end
end

% --------------------- 输出最终结果 ---------------------
fprintf('\n=== 对比完成 ===\n');
for si = 1:length(systems)
    fprintf('%s RMS误差：%.3f 米\n', sys_names{si}, all_rms(si));
end