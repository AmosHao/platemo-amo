%% ========== 直接运行区域：修改下面的参数后直接运行脚本 ==========
% 直接运行时的默认参数（修改这里即可）
excel_file = '20260214order.xlsx';  % Excel文件名
obj_col = 2;  % obj列（1=f1总能量, 2=f2总时间）
rank_idx = 1;  % 索引（第几小的值，如1表示最小值）

% 调用函数绘制
plot_order_map_path_func(excel_file, obj_col, rank_idx);

%% ========== 函数定义 ==========
function plot_order_map_path_func(excel_file, obj_col, rank_idx)
% 绘制订单配送问题地图并显示指定解的路径
% 
% 参数:
%   excel_file - Excel文件名（在 Order_Data/testData 下）或完整路径
%   obj_col    - 目标列（1=f1总能量, 2=f2总时间）
%   rank_idx   - 索引（第几小的值，如1表示最小值）
%
% 示例:
%   plot_order_map_path_func('20260214order_obj3minf1.xlsx', 1, 1);  % 绘制f1最小的解的路径
%   plot_order_map_path_func('20260214order_obj3minf1.xlsx', 2, 1);  % 绘制f2最小的解的路径

    % 默认参数
    if nargin < 1 || isempty(excel_file)
        excel_file = '20260214order_obj3minf1.xlsx';
    end
    if nargin < 2 || isempty(obj_col)
        obj_col = 1;  % 默认f1列
    end
    if nargin < 3 || isempty(rank_idx)
        rank_idx = 1;  % 默认最小值
    end
    
    % 引入数据文件
    order_data;
    
    % 读取Excel文件
    script_dir = fileparts(mfilename('fullpath'));
    [filepath, ~, ~] = fileparts(excel_file);
    if isempty(filepath)
        excel_file = fullfile(script_dir, '..', 'Order_Data', 'testData', excel_file);
    end
    
    if ~isfile(excel_file)
        error('文件不存在: %s', excel_file);
    end
    
    % 读取obj_data和dec_data
    obj_data = readmatrix(excel_file, 'Sheet', 'obj_data');
    dec_data = readmatrix(excel_file, 'Sheet', 'dec_data');
    
    if isempty(obj_data) || isempty(dec_data)
        error('无法读取Excel数据，请检查Sheet名称');
    end
    
    % 找到指定obj列的第rank_idx小的值对应的索引
    obj_col_data = obj_data(:, obj_col);
    valid_idx = isfinite(obj_col_data);
    obj_col_data_valid = obj_col_data(valid_idx);
    [sorted_obj, sort_idx] = sort(obj_col_data_valid);
    
    if rank_idx > length(sorted_obj)
        error('索引 %d 超出范围，有效解数量为 %d', rank_idx, length(sort_idx));
    end
    
    % 找到原始索引（考虑NaN行）
    valid_positions = find(valid_idx);
    target_idx = valid_positions(sort_idx(rank_idx));
    
    % 获取对应的dec值
    if target_idx > size(dec_data, 1)
        error('索引 %d 超出dec_data范围', target_idx);
    end
    dec_sequence = dec_data(target_idx, :);
    
    % 显示信息
    fprintf('选择解: obj列%d的第%d小值 = %.2f, 索引 = %d\n', ...
        obj_col, rank_idx, sorted_obj(rank_idx), target_idx);

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

% 设置坐标轴范围（路径点在 dotss 中 z=800m，需包含在内才能看到路径）
xlim([0, 10000]);
ylim([0, 10000]);
zlim([0, 1000]);  % 0-200m 障碍物，800m 路径点与商家/客户/配送中心，故取 0-1000

grid on;

% ========== 绘制无人机路径（如果提供了dec序列） ==========
if exist('dec_sequence', 'var') && ~isempty(dec_sequence)
    % 解析dec序列：用0分隔成多段
    zero_positions = find(dec_sequence == 0);
    segments = {};
    
    if isempty(zero_positions)
        % 没有0分隔符，整个序列是一段
        segments{1} = dec_sequence(dec_sequence ~= 0);
    else
        % 用0分隔成多段
        start_idx = 1;
        for i = 1:length(zero_positions)
            end_idx = zero_positions(i) - 1;
            if end_idx >= start_idx
                seg = dec_sequence(start_idx:end_idx);
                seg = seg(seg ~= 0);  % 移除可能的0
                if ~isempty(seg)
                    segments{end+1} = seg;
                end
            end
            start_idx = zero_positions(i) + 1;
        end
        % 最后一段（0之后）
        if start_idx <= length(dec_sequence)
            seg = dec_sequence(start_idx:end);
            seg = seg(seg ~= 0);
            if ~isempty(seg)
                segments{end+1} = seg;
            end
        end
    end
    
    % 为每架无人机分配不同颜色
    colors = lines(length(segments));  % 使用lines颜色映射
    
    % 绘制每架无人机的路径
    for seg_idx = 1:length(segments)
        seg = segments{seg_idx};
        if isempty(seg)
            continue;
        end
        
        % 构建完整路径：配送中心 -> 段内点 -> 配送中心
        path_indices = [21, seg, 21];  % 21是配送中心索引
        
        % 获取路径点的坐标
        path_coords = zeros(length(path_indices), 3);
        for i = 1:length(path_indices)
            idx = path_indices(i);
            if idx >= 1 && idx <= size(dotss, 1)
                path_coords(i, :) = dotss(idx, :);
            else
                warning('索引 %d 超出范围，跳过', idx);
                path_coords(i, :) = [NaN, NaN, NaN];
            end
        end
        
        % 移除NaN点
        valid_path = ~isnan(path_coords(:, 1));
        path_coords = path_coords(valid_path, :);
        
        if size(path_coords, 1) < 2
            continue;
        end
        
        % 绘制路径连线（使用不同颜色，使用dotss中的z坐标）
        plot3(path_coords(:, 1), path_coords(:, 2), path_coords(:, 3), ...
            'Color', colors(seg_idx, :), 'LineWidth', 2.5, 'LineStyle', '-', ...
            'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(seg_idx, :));
        
        % 添加无人机标签
        if ~isempty(path_coords)
            mid_point = path_coords(ceil(size(path_coords, 1)/2), :);
            text(mid_point(1), mid_point(2), mid_point(3) + 100, ...
                sprintf('无人机%d', seg_idx), 'FontSize', 9, ...
                'Color', colors(seg_idx, :), 'FontWeight', 'bold');
        end
    end
    
    % 更新图例，添加无人机路径
    for seg_idx = 1:length(segments)
        h_path = plot3(NaN, NaN, NaN, 'Color', colors(seg_idx, :), ...
            'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 8, ...
            'MarkerFaceColor', colors(seg_idx, :));
        legend_handles = [legend_handles, h_path];
        legend_items{end+1} = sprintf('无人机%d路径', seg_idx);
    end
    
    % 更新图例
    legend(legend_handles, legend_items, 'Location', 'best');
    
    % 更新标题
    obj_name = {'总能量', '总时间'};
    title(sprintf('订单配送问题地图 - %s第%d小解 (索引%d)', ...
        obj_name{obj_col}, rank_idx, target_idx));
end

hold off;
end
