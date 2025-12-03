%% --------------------- 地面投影对比（最终修正版） ---------------------
figure('Position',[100 100 1000 500], 'Color','w');
hold on; axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)');
title('Ground Track Comparison (2019-12-01)');

% 根据图片实际显示修正颜色
% 1. GPS MEO（紫色虚线） - 图片中实际是紫色虚线
gps_color = [0.5 0 0.5]; % 紫色
h_gps = plot([], [], 'Color', gps_color, 'LineStyle', '--', 'LineWidth', 1);
for i = 1:numel(prns_gps)
    idx = gps.PRN == prns_gps(i);
    plot(gps.X(idx), gps.Y(idx), 'Color', gps_color, 'LineStyle', '--', 'LineWidth', 1);
end

% 2. BDS GEO（蓝色实线） - 图片中实际是蓝色实线
bds_geo_color = [0 0 0.8]; % 深蓝色
h_geo = plot([], [], 'Color', bds_geo_color, 'LineStyle', '-', 'LineWidth', 3);
for prn = 1:5
    idx = bds.PRN == prn;
    if any(idx)
        plot(bds.X(idx), bds.Y(idx), 'Color', bds_geo_color, 'LineStyle', '-', 'LineWidth', 3);
    end
end

% 3. BDS IGSO（绿色实线） - 图片中实际是绿色实线
bds_igso_color = 'g'; % 绿色
h_igso = plot([], [], 'Color', bds_igso_color, 'LineStyle', '-', 'LineWidth', 1.8);
for prn = 6:14
    idx = bds.PRN == prn;
    if any(idx)
        plot(bds.X(idx), bds.Y(idx), 'Color', bds_igso_color, 'LineStyle', '-', 'LineWidth', 1.8);
    end
end

% 4. BDS MEO（深蓝色实线） - 图片中实际是深蓝色实线
bds_meo_color = [0 0 0.5]; % 更深的蓝色
h_meo = plot([], [], 'Color', bds_meo_color, 'LineStyle', '-', 'LineWidth', 1.2);
for prn = 15:30
    idx = bds.PRN == prn;
    if any(idx)
        plot(bds.X(idx), bds.Y(idx), 'Color', bds_meo_color, 'LineStyle', '-', 'LineWidth', 1.2);
    end
end

% 修正图例显示
legend([h_gps, h_geo, h_igso, h_meo], ...
       {'GPS MEO', 'BDS GEO', 'BDS IGSO', 'BDS MEO'}, ...
       'Location', 'bestoutside', ...
       'Color', 'w', 'EdgeColor', 'k');
hold off;