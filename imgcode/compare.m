%% ==================== 2. 广播星历 vs 精密星历对比 ====================
close all; clear; clc;

% ---------- 1. 读取你自己算出来的广播星历文件（15 min） ----------
% 文件格式示例：
% * 2019 12 1 0 15 0.00000000 900.000
% G01 32715679.123 1234568.123 5432168.123 0.00012456
% ...

brdc_gps = readBrdcSP3('gps_position.txt');   % GPS 广播
brdc_bds = readBrdcSP3('gps_position.txt');   % 北斗广播

% ---------- 2. 读取 IGS/COD 精密星历（.sp3） ----------
% 下载地址示例（2019-12-01 为 GPS 周 2082 第0天）：
% https://cddis.nasa.gov/archive/gnss/products/2082/igs20820.sp3.Z
% 解压后得到 igs20820.sp3
prec = readSP3('WUM0MGXFIN_20193350000_01D_15M_ORB.SP3');   % 同时包含 GPS + BDS（COD/IGS 都带）

% ---------- 3. 对齐历元并计算误差 ----------
systems = {'G','C'};   % GPS 和 BDS
all_rms = [];

figure('Position',[100 100 1100 700]);
tiledlayout(2,1,'TileSpacing','compact');

for si = 1:2
    sys = systems{si};
    if sys=='G', brdc = brdc_gps; else, brdc = brdc_bds; end
    
    prns = unique(brdc.PRN);
    err_all = [];   % 保存所有误差（m）
    epoch_idx = 0;
    
    for i = 1:height(brdc)
        epoch_idx = epoch_idx + 1;
        prn = brdc.PRN(i);
        Xb = [brdc.X(i), brdc.Y(i), brdc.Z(i)] * 1000;   % 广播星历是 km → m
        
        % 在精密星历里找同一颗星、同一历元（时间误差<1s就算匹配）
        idx_prec = find(prec.PRN==prn & ...
                        abs(prec.t - brdc.t(i)) < 1/86400);   % 1秒以内
        
        if ~isempty(idx_prec)
            Xp = [prec.X(idx_prec), prec.Y(idx_prec), prec.Z(idx_prec)] * 1000; % km→m
            err = norm(Xb - Xp);           % 3D 轨道误差（m）
            err_all(epoch_idx, prn) = err;
        else
            err_all(epoch_idx, prn) = NaN;
        end
    end
    
    % ----------- 画每颗星的误差时序图 -----------
    nexttile;
    hold on; grid on;
    colors = lines(numel(prns));
    for k = 1:numel(prns)
        prn = prns(k);
        plot(err_all(:,prn), 'Color',colors(k,:), 'LineWidth',1.2, ...
             'DisplayName',sprintf('%s%02d',sys,prn));
    end
    xlabel('Epoch (15 min interval)');
    ylabel('Orbit Error (m)');
    title(sprintf('%s Broadcast vs Precise Ephemeris Error', ...
                  ifelse(sys=='G','GPS','BDS')));
    legend('show','Location','eastoutside');
    
    % ----------- 计算并显示 RMS -----------
    valid = ~isnan(err_all(:));
    rms_total = sqrt(mean(err_all(valid).^2));
    fprintf('%s 总体轨道 RMS = %.3f m\n', ifelse(sys=='G','GPS','BDS'), rms_total);
    all_rms = [all_rms; rms_total];
end

disp('=== 全部完成 ===');

function data = readBrdcSP3(filename)
% 读取你自己生成的类似SP3格式的广播星历文件
fid = fopen(filename,'r');
if fid==-1, error('文件打不开: %s',filename); end

t = []; prn = []; X = []; Y = []; Z = []; clk = [];

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    
    if startsWith(line,'*')   % 历元行
        parts = split(line);
        year = str2double(parts{2}); month = str2double(parts{3});
        day  = str2double(parts{4}); hour = str2double(parts{5});
        minute = str2double(parts{6}); sec = str2double(parts{7});
        % 转成 MATLAB datenum（便于后面比较）
        epoch_time = datenum([year month day hour minute sec]);
    elseif startsWith(line,'G') || startsWith(line,'C')
        sys = line(1);
        prn_num = str2double(line(2:3));
        vals = sscanf(line(5:end),'%f%f%f%f');   % X Y Z clk(秒)
        t(end+1,1)    = epoch_time;
        prn(end+1,1)  = prn_num + ifelse(sys=='C',100,0); % C用100+prn 区分
        X(end+1,1)    = vals(1)/1000;   % 原始是 m → km（后面统一转）
        Y(end+1,1)    = vals(2)/1000;
        Z(end+1,1)    = vals(3)/1000;
        clk(end+1,1)  = vals(4);
    end
end
fclose(fid);

data = table(t, prn, X, Y, Z, clk, 'VariableNames',{'t','PRN','X','Y','Z','clk'});
end

function data = readSP3(filename)
% 读取标准 IGS/COD .sp3 文件（支持 GPS + BDS）
fid = fopen(filename,'r');
if fid==-1, error('打不开精密星历: %s',filename); end

t = []; prn = []; X = []; Y = []; Z = []; clk = [];

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    
    if startsWith(line,'*')   % 历元行
        parts = split(line);
        year = str2double(parts{2}); month = str2double(parts{3});
        day  = str2double(parts{4}); hour = str2double(parts{5});
        minute = str2double(parts{6}); sec = str2double(parts{7});
        epoch_time = datenum([year month day hour minute sec]);
    elseif startsWith(line,'P') && (line(2)=='G' || line(2)=='C')
        sys = line(2);
        prn_num = str2double(line(3:4));
        vals = sscanf(line(5:end),'%f%f%f%f');
        if numel(vals)>=4
            t(end+1,1) = epoch_time;
            prn(end+1,1) = prn_num + ifelse(sys=='C',100,0);
            X(end+1,1) = vals(1);      % km
            Y(end+1,1) = vals(2);
            Z(end+1,1) = vals(3);
            clk(end+1,1) = vals(4)*1e-6;   % 通常是微秒 → 秒
        end
    end
end
fclose(fid);

data = table(t, prn, X, Y, Z, clk, ...
            'VariableNames',{'t','PRN','X','Y','Z','clk'});
end

function out = ifelse(cond,a,b)
if cond, out=a; else, out=b; end
end