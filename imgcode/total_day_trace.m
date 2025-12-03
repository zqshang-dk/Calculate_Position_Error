%% ==================== 1. 分别绘制 GPS 和 BDS 全天轨迹（两张图）===================
close all; clear; clc;

% --------------------- 读取每60s的轨迹CSV ---------------------
gps = readtable('gps_track.csv');   % 列: PRN, t, X, Y, Z
bds = readtable('bds_track.csv');   % 同上

Rearth = 6378137;  % WGS84 平均半径（m）

%% --------------------- 图1：GPS 卫星轨迹 ---------------------
figure('Position', [100, 100, 900, 800], 'Color', 'white');
hold on; grid on; axis equal;
xlabel('X (m)', 'FontSize', 12); 
ylabel('Y (m)', 'FontSize', 12); 
zlabel('Z (m)', 'FontSize', 12);
title('GPS Satellite Orbits (Broadcast Ephemeris, 2019-12-01, 60s interval)', ...
      'FontSize', 14, 'FontWeight', 'bold');

% 画地球
[Xs,Ys,Zs] = sphere(50);
h_earth = surf(Rearth*Xs, Rearth*Ys, Rearth*Zs, ...
    'FaceColor', [0.8 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
camlight headlight; lighting gouraud;

% GPS 轨迹
prns_gps = unique(gps.PRN);
colors = lines(numel(prns_gps));

for i = 1:numel(prns_gps)
    idx = gps.PRN == prns_gps(i);
    plot3(gps.X(idx), gps.Y(idx), gps.Z(idx), ...
        'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('G%02d', prns_gps(i)));
end

legend('show', 'Location', 'bestoutside', 'FontSize', 10);
view(45, 25);
xlim([-3e7, 3e7]); ylim([-3e7, 3e7]); zlim([-3e7, 3e7]);
set(gca, 'FontSize', 11);
hold off;

%% --------------------- 图2：BDS 卫星轨迹 ---------------------
figure('Position', [100, 100, 900, 800], 'Color', 'white');
hold on; grid on; axis equal;
xlabel('X (m)', 'FontSize', 12); 
ylabel('Y (m)', 'FontSize', 12); 
zlabel('Z (m)', 'FontSize', 12);
title('BDS Satellite Orbits (Broadcast Ephemeris, 2019-12-01, 60s interval)', ...
      'FontSize', 14, 'FontWeight', 'bold');

% 画地球
h_earth = surf(Rearth*Xs, Rearth*Ys, Rearth*Zs, ...
    'FaceColor', [0.8 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
camlight headlight; lighting gouraud;

% BDS 轨迹（不同类型用不同线型）
prns_bds = unique(bds.PRN);
colors = lines(numel(prns_bds));

for i = 1:numel(prns_bds)
    prn = prns_bds(i);
    idx = bds.PRN == prn;
    
    % 判断卫星类型（粗略：PRN <= 5 为 GEO，其余为 MEO/IGSO）
    if prn <= 5
        % GEO：红色实线，较粗
        plot3(bds.X(idx), bds.Y(idx), bds.Z(idx), ...
            'r-', 'LineWidth', 2.5, 'DisplayName', sprintf('C%02d (GEO)', prn));
    elseif prn <= 18
        % IGSO：蓝色虚线
        plot3(bds.X(idx), bds.Y(idx), bds.Z(idx), ...
            'b--', 'LineWidth', 1.8, 'DisplayName', sprintf('C%02d (IGSO)', prn));
    else
        % MEO：绿色实线
        plot3(bds.X(idx), bds.Y(idx), bds.Z(idx), ...
            'g-', 'LineWidth', 1.5, 'DisplayName', sprintf('C%02d (MEO)', prn));
    end
end

legend('show', 'Location', 'bestoutside', 'FontSize', 10);
view(120, 30);  % 从东经105°附近看，更能看出 GEO 分布
xlim([-5e7, 5e7]); ylim([-5e7, 5e7]); zlim([-5e7, 5e7]);
set(gca, 'FontSize', 11);
hold off;

%% --------------------- 可选：保存高清图 ---------------------
% print(1, 'GPS_Orbit_3D.png', '-dpng', '-r600');
% print(2, 'BDS_Orbit_3D.png', '-dpng', '-r600');