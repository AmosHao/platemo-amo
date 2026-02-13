% 绘制订单配送问题地图
% 包含：配送中心、商家、客户点、禁飞区、人流密集区、噪音敏感区

% 引入数据文件
order_data;

% 创建图形窗口
figure;
hold on;
view(3);
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('订单配送问题地图');

% ========== 绘制禁飞区（矩形，红色） ==========
for i = 1:size(forbidden_zones, 1)
    rect = forbidden_zones(i, :);
    x_min = rect(1);
    y_min = rect(2);
    x_max = rect(3);
    y_max = rect(4);
    zValue = 0;
    
    % 绘制矩形底面
    fill3([x_min, x_max, x_max, x_min], ...
          [y_min, y_min, y_max, y_max], ...
          [zValue, zValue, zValue, zValue], ...
          [1, 0, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2);  % 红色，半透明
    
    % 绘制矩形顶面（高度设为120米）
    z_top = 120;
    fill3([x_min, x_max, x_max, x_min], ...
          [y_min, y_min, y_max, y_max], ...
          [z_top, z_top, z_top, z_top], ...
          [1, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 1);
    
    % 绘制矩形侧面
    % 前面
    fill3([x_min, x_max, x_max, x_min], ...
          [y_min, y_min, y_min, y_min], ...
          [zValue, zValue, z_top, z_top], ...
          [1, 0, 0], 'FaceAlpha', 0.3);
    % 后面
    fill3([x_min, x_max, x_max, x_min], ...
          [y_max, y_max, y_max, y_max], ...
          [zValue, zValue, z_top, z_top], ...
          [1, 0, 0], 'FaceAlpha', 0.3);
    % 左面
    fill3([x_min, x_min, x_min, x_min], ...
          [y_min, y_max, y_max, y_min], ...
          [zValue, zValue, z_top, z_top], ...
          [1, 0, 0], 'FaceAlpha', 0.3);
    % 右面
    fill3([x_max, x_max, x_max, x_max], ...
          [y_min, y_max, y_max, y_min], ...
          [zValue, zValue, z_top, z_top], ...
          [1, 0, 0], 'FaceAlpha', 0.3);
end

% ========== 绘制人流密集区（圆形，淡绿色） ==========
for i = 1:size(crowded_zones, 1)
    circle = crowded_zones(i, :);
    cx = circle(1);
    cy = circle(2);
    radius = circle(3);
    
    % 生成圆形数据
    theta = linspace(0, 2*pi, 100);
    xCircle = cx + radius * cos(theta);
    yCircle = cy + radius * sin(theta);
    zCircle = zeros(size(theta));
    
    % 绘制圆形底面
    fill3(xCircle, yCircle, zCircle, [0.6, 1, 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'LineWidth', 1.5);
    
    % 绘制圆形顶面（高度设为120米）
    z_top = 120;
    fill3(xCircle, yCircle, z_top * ones(size(theta)), [0.6, 1, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'g', 'LineWidth', 1);
    
    % 绘制圆柱侧面
    [X, Y, Z] = cylinder(radius, 100);
    X = X + cx;
    Y = Y + cy;
    Z = Z * z_top;
    surf(X, Y, Z, 'FaceColor', [0.6, 1, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% ========== 绘制噪音敏感区（矩形，橙色） ==========
for i = 1:size(noise_zones, 1)
    rect = noise_zones(i, :);
    x_min = rect(1);
    y_min = rect(2);
    x_max = rect(3);
    y_max = rect(4);
    zValue = 0;
    
    % 绘制矩形底面
    fill3([x_min, x_max, x_max, x_min], ...
          [y_min, y_min, y_max, y_max], ...
          [zValue, zValue, zValue, zValue], ...
          [1, 0.5, 0], 'FaceAlpha', 0.4, 'EdgeColor', [1, 0.5, 0], 'LineWidth', 1.5);  % 橙色，半透明
    
    % 绘制矩形顶面（高度设为120米）
    z_top = 120;
    fill3([x_min, x_max, x_max, x_min], ...
          [y_min, y_min, y_max, y_max], ...
          [z_top, z_top, z_top, z_top], ...
          [1, 0.5, 0], 'FaceAlpha', 0.2, 'EdgeColor', [1, 0.5, 0], 'LineWidth', 1);
    
    % 绘制矩形侧面
    fill3([x_min, x_max, x_max, x_min], ...
          [y_min, y_min, y_min, y_min], ...
          [zValue, zValue, z_top, z_top], ...
          [1, 0.5, 0], 'FaceAlpha', 0.2);
    fill3([x_min, x_max, x_max, x_min], ...
          [y_max, y_max, y_max, y_max], ...
          [zValue, zValue, z_top, z_top], ...
          [1, 0.5, 0], 'FaceAlpha', 0.2);
    fill3([x_min, x_min, x_min, x_min], ...
          [y_min, y_max, y_max, y_min], ...
          [zValue, zValue, z_top, z_top], ...
          [1, 0.5, 0], 'FaceAlpha', 0.2);
    fill3([x_max, x_max, x_max, x_max], ...
          [y_min, y_max, y_max, y_min], ...
          [zValue, zValue, z_top, z_top], ...
          [1, 0.5, 0], 'FaceAlpha', 0.2);
end

% ========== 绘制配送中心（紫色大圆点） ==========
depot = dotss(21, :);  % 配送中心在第21行（最后一行）
% 配送中心高度设为0（地面）
scatter3(depot(1), depot(2), 0, 200, [0.8, 0.6, 0.8], 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 2);
% 标签位置：在点的上方（y方向），z=250，避免遮挡
text(depot(1), depot(2) + 200, 250, '配送中心', 'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% ========== 绘制商家（蓝色圆点，编号1-10） ==========
merchants = dotss(1:10, :);  % 商家在第1-10行
% 商家高度设为0（地面）
scatter3(merchants(:, 1), merchants(:, 2), zeros(size(merchants, 1), 1), 100, 'b', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 添加商家标签（标签位置在点的上方y方向，z=250，避免遮挡）
for i = 1:size(merchants, 1)
    text(merchants(i, 1), merchants(i, 2) + 200, 250, ...
        sprintf('商家%d', i), 'FontSize', 8, 'Color', 'b', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% ========== 绘制客户点（绿色圆点，编号11-20） ==========
customers = dotss(11:20, :);  % 客户点在第11-20行
% 客户点高度设为0（地面）
scatter3(customers(:, 1), customers(:, 2), zeros(size(customers, 1), 1), 100, 'g', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 添加客户点标签（标签位置在点的上方y方向，z=250，避免遮挡）
for i = 1:size(customers, 1)
    text(customers(i, 1), customers(i, 2) + 200, 250, ...
        sprintf('客户%d', i+10), 'FontSize', 8, 'Color', 'g', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% ========== 添加图例 ==========
legend_items = {};
legend_handles = [];

% 配送中心
h1 = scatter3(NaN, NaN, NaN, 200, [0.8, 0.6, 0.8], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
legend_handles = [legend_handles, h1];
legend_items{end+1} = '配送中心';

% 商家
h2 = scatter3(NaN, NaN, NaN, 100, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_handles = [legend_handles, h2];
legend_items{end+1} = '商家';

% 客户点
h3 = scatter3(NaN, NaN, NaN, 100, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_handles = [legend_handles, h3];
legend_items{end+1} = '客户点';

% 禁飞区
h4 = fill3(NaN, NaN, NaN, [1, 0, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2);
legend_handles = [legend_handles, h4];
legend_items{end+1} = '禁飞区';

% 人流密集区
h5 = fill3(NaN, NaN, NaN, [0.6, 1, 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'LineWidth', 1.5);
legend_handles = [legend_handles, h5];
legend_items{end+1} = '人流密集区';

% 噪音敏感区
h6 = fill3(NaN, NaN, NaN, [1, 0.5, 0], 'FaceAlpha', 0.4, 'EdgeColor', [1, 0.5, 0], 'LineWidth', 1.5);
legend_handles = [legend_handles, h6];
legend_items{end+1} = '噪音敏感区';

legend(legend_handles, legend_items, 'Location', 'best');

% 设置坐标轴范围
xlim([0, 10000]);
ylim([0, 10000]);
zlim([0, 200]);  % z轴范围设为0-200m，确保120m的障碍物和150m的标签都能显示

grid on;
hold off;

