%% ========== 直接运行区域：修改下面的参数后直接运行脚本 ==========
% 两种模式：
% - Excel 模式：直接给 excel_file（老逻辑）
% - Ablation 模式：给 ablation_dir + expName + hv_rank（先按 HV 排名选 run，再在该 run 内按目标值排名选解）
%
% Ablation 数据结构与 runplatemo_amo_EXP_m.m 一致（仅写 .mat，不写单 run 的 xlsx）：
%   <ablation_dir>/<expName>/<expName>_runs.mat
%   - runs(r).Run, .hv_end, .Obj_end, .Dec_end, .Details, ...
%   - summary（可选，绘图以 runs 为准）

mode = "ablation";   % "ablation" | "excel"

% --- Ablation 模式参数 ---
ablation_dir = 'D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\newdata2604\ablation0416';
expName      = 'ablation0416_noJF_NSGAIII_G';  % 例如：ablation0415_NSGAII_G
hv_rank      = 1;    % 在该实验的 31 次（或实际次数）里按 HV_end 降序的名次：1=最好

% --- 选解参数（两种模式通用）---
obj_col         = 1;    % 目标列（1=f1总能量, 2=f2总时间）
rank_idx        = 1;    % 该目标列的第几小（1=最小）
filter_feasible = true; % 是否只从可行解中选取（true=筛出details_data里 Feasible=1 的解）

% --- Excel 模式参数 ---
excel_file = '2026022224_09order.xlsx';  % Excel文件名（在 Order_Data/testData 下）或完整路径

% --- 仅地图（无路径）---
% true：只画禁飞区/人流区/噪声区、配送中心、商家与客户点，并显示图例；不读消融或 Excel 解、不画路径、无图题
map_only_env = true;

% 调用入口
plot_order_map_path_select_and_draw(mode, excel_file, ablation_dir, expName, hv_rank, obj_col, rank_idx, filter_feasible, map_only_env);

%% ========== 函数定义 ==========
function plot_order_map_path_select_and_draw(mode, excel_file, ablation_dir, expName, hv_rank, obj_col, rank_idx, filter_feasible, map_only_env)
% 绘制订单配送问题地图并显示指定解的路径
% 
% 参数:
%   mode            - "excel" 或 "ablation"
%   excel_file      - Excel文件名（在 Order_Data/testData 下）或完整路径（excel 模式用）
%   ablation_dir    - 消融实验根目录（ablation 模式用）
%   expName         - 实验子文件夹名（ablation 模式用，例如 ablation0415_NSGAII_G）
%   hv_rank         - HV_end 排名名次（ablation 模式用，1=最好）
%   obj_col         - 目标列（1=f1总能量, 2=f2总时间）
%   rank_idx        - 索引（第几小的值，如1表示最小值）
%   filter_feasible - 是否只从可行解中选取（true 时读 Details/details_data
%                     最后一列 Feasible=1，再按 obj_col/rank_idx 排序选解）
%   map_only_env    - true：仅环境与节点图+图例，无路径、无标题、不加载解
%
% 示例:
%   plot_order_map_path_select_and_draw("excel",'20260222order.xlsx',[],[],[], 1,1,true,false);
%   plot_order_map_path_select_and_draw("ablation",[], '...\ablation0415','ablation0415_NSGAII_G', 1, 2, 1, true, false);
%   plot_order_map_path_select_and_draw("ablation",[], '...\ablation0415','ablation0415_NSGAII_G', 1, 1, 1, false, true);

    % 默认参数
    if nargin < 1 || isempty(mode)
        mode = "excel";
    end
    if nargin < 2
        excel_file = '';
    end
    if nargin < 3
        ablation_dir = '';
    end
    if nargin < 4
        expName = '';
    end
    if nargin < 5 || isempty(hv_rank)
        hv_rank = 1;
    end
    if nargin < 6 || isempty(obj_col)
        obj_col = 1;  % 默认f1列
    end
    if nargin < 7 || isempty(rank_idx)
        rank_idx = 1;  % 默认最小值
    end
    if nargin < 8 || isempty(filter_feasible)
        filter_feasible = false;
    end
    if nargin < 9 || isempty(map_only_env)
        map_only_env = false;
    end
    map_only_env = logical(map_only_env);
    
    % 引入数据文件
    order_data;

    dec_sequence = [];  % map_only_env 为 true 时不读取解，保持为空则不画路径
    if ~map_only_env
    % 根据模式得到 obj_data / dec_data
    % Ablation（runplatemo_amo_EXP_m）：权威来源为 <expName>_runs.mat 内的 runs(r)
    mode = string(mode);
    details_mat = [];  % 非空时表示可行解信息来自 mat 的 Details 矩阵
    if mode == "ablation"
        expDir = fullfile(ablation_dir, expName);
        if ~isfolder(expDir)
            error('实验目录不存在: %s', expDir);
        end
        summaryPath = fullfile(expDir, [expName '_summary.xlsx']);
        matPath     = fullfile(expDir, [expName '_runs.mat']);
        if isfile(matPath)
            % 与 runplatemo_amo_EXP_m 一致：按 runs(:).hv_end 降序排名
            [runIdx, hv_end] = get_run_index_by_HV_rank_mat(matPath, hv_rank);
            fprintf('[Ablation] 选择实验 %s：HV 排名第 %d -> run%d (HV_end=%.6g) [来源: %s]\n', ...
                expName, hv_rank, runIdx, hv_end, matPath);
            excel_file = matPath;
            [obj_data, dec_data, details_mat] = load_run_from_runs_mat(matPath, runIdx);
        elseif isfile(summaryPath)
            % 旧流程兜底：仅有 summary xlsx 时再尝试按表排名 + 找单 run xlsx
            warning('myapp:noRunsMat', ...
                '未找到 %s（runplatemo_amo_EXP_m 标准输出）。将尝试 _summary.xlsx + run xlsx。', matPath);
            [runIdx, hv_end] = get_run_index_by_HV_rank_xlsx(summaryPath, hv_rank);
            fprintf('[Ablation] 选择实验 %s：HV 排名第 %d -> run%d (HV_end=%.6g) [来源: summary xlsx]\n', ...
                expName, hv_rank, runIdx, hv_end);
            runXlsx = find_run_excel_file(expDir, expName, runIdx);
            if isempty(runXlsx)
                error(['缺少 %s，且目录下无对应 run 的 xlsx。\n' ...
                    '请确认已运行 runplatemo_amo_EXP_m 并生成 *_runs.mat。'], matPath);
            end
            excel_file = runXlsx;
            obj_data = readmatrix(excel_file, 'Sheet', 'obj_data');
            dec_data = readmatrix(excel_file, 'Sheet', 'dec_data');
        else
            error(['未找到消融结果文件。\n' ...
                'runplatemo_amo_EXP_m 应生成: %s\n' ...
                '（可选旧格式）: %s'], matPath, summaryPath);
        end
    else
        % Excel 模式：按原逻辑补全相对路径
        if isempty(excel_file)
            excel_file = '20260214order_obj3minf1.xlsx';
        end
        script_dir = fileparts(mfilename('fullpath'));
        [filepath, ~, ~] = fileparts(excel_file);
        if isempty(filepath)
            excel_file = fullfile(script_dir, '..', 'Order_Data', 'testData', excel_file);
        end
        if ~isfile(excel_file)
            error('文件不存在: %s', excel_file);
        end
        obj_data = readmatrix(excel_file, 'Sheet', 'obj_data');
        dec_data = readmatrix(excel_file, 'Sheet', 'dec_data');
    end
    
    if isempty(obj_data) || isempty(dec_data)
        error('无法读取目标数据：obj_data/dec_data 为空（来源: %s）', excel_file);
    end
    
    % ---- 筛选可行解 ----
    % 读取 details_data，提取 Feasible=1（最后一列）的行索引，
    % 或将 obj_data / dec_data 限定在这些行上，再做排序选解
    feasible_mask = true(size(obj_data, 1), 1);  % 默认全部保留
    if filter_feasible
        got_details = false;
        if ~isempty(details_mat)
            details = details_mat;
            got_details = true;
        else
            try
                details = readmatrix(excel_file, 'Sheet', 'details_data');
                got_details = ~isempty(details);
            catch e
                warning('myapp:readDetailsFailed', '读取 details_data 失败（%s），忽略可行解筛选', e.message);
            end
        end
        if got_details
            feasible_col = details(:, end);   % 最后一列为 Feasible 标志
            n_rows = min(size(obj_data, 1), numel(feasible_col));
            feasible_mask = false(size(obj_data, 1), 1);
            feasible_mask(1:n_rows) = (feasible_col(1:n_rows) == 1);
            n_feasible = sum(feasible_mask);
            fprintf('[筛选] details 共 %d 行，可行解 %d 个\n', size(details,1), n_feasible);
            if n_feasible == 0
                warning('myapp:noFeasible', '%s', '没有找到任何可行解（Feasible=1），将使用全部解');
                feasible_mask = true(size(obj_data, 1), 1);
            end
        else
            warning('myapp:emptyDetails', '%s', '无 details 数据，忽略可行解筛选');
        end
    end
    
    % 找到指定obj列的第rank_idx小的值对应的索引
    % 先在可行掩码范围内挑选有效行（finite + feasible）
    obj_col_data = obj_data(:, obj_col);
    valid_idx = isfinite(obj_col_data) & feasible_mask;
    obj_col_data_valid = obj_col_data(valid_idx);
    [sorted_obj, sort_idx] = sort(obj_col_data_valid);
    
    if rank_idx > length(sorted_obj)
        error('索引 %d 超出范围，有效解数量为 %d', rank_idx, length(sort_idx));
    end
    
    % 找到原始索引（考虑NaN行及feasible掩码）
    valid_positions = find(valid_idx);
    target_idx = valid_positions(sort_idx(rank_idx));
    
    % 获取对应的dec值
    if target_idx > size(dec_data, 1)
        error('索引 %d 超出dec_data范围', target_idx);
    end
    dec_sequence = dec_data(target_idx, :);
    
    % 显示信息
    feasible_tag = '';
    if filter_feasible
        feasible_tag = '（仅可行解）';
    end
    % 打印目标函数值（同时打印 f1/f2，便于核对）
    obj_row = obj_data(target_idx, :);
    f1 = NaN; f2 = NaN;
    if numel(obj_row) >= 1, f1 = obj_row(1); end
    if numel(obj_row) >= 2, f2 = obj_row(2); end
    fprintf('选择解%s: obj列%d的第%d小值 = %.6g, 索引 = %d\n', ...
        feasible_tag, obj_col, rank_idx, sorted_obj(rank_idx), target_idx);
    fprintf('  目标函数值: f1=%.6g, f2=%.6g', f1, f2);
    if obj_col <= numel(obj_row)
        fprintf(' （当前排序列 obj_col=%d 的值=%.6g）\n', obj_col, obj_row(obj_col));
    else
        fprintf('\n');
    end
    else
        fprintf('[map_only_env] 仅绘制环境、商家/客户点与图例；不加载解数据、不画路径、无图题。\n');
    end

% 创建图形窗口（平面地图：坐标以 km 显示；尺寸适中，map_only 时通过轴位置为右侧图例留边）
figH = figure('Color', 'w', 'Units', 'pixels', 'Position', [80 80 920 720]);
hold on;
m2km = 1e-3;
forbidden_zones_km = forbidden_zones * m2km;
noise_zones_km     = noise_zones * m2km;
crowded_zones_km   = crowded_zones;
crowded_zones_km(:, 1:3) = crowded_zones(:, 1:3) * m2km;
dotss_km = [dotss(:, 1) * m2km, dotss(:, 2) * m2km, dotss(:, 3) * m2km];
z_top_obstacle = 120 * m2km;

view(2);
axis equal;
xlabel('X/km');
ylabel('Y/km');
% 图题在全部绘制完成后设置（仅实验名 expName），见 hold off 前

% ========== 绘制禁飞区（矩形，红色） ==========
for i = 1:size(forbidden_zones_km, 1)
    rect = forbidden_zones_km(i, :);
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
    
    % 绘制矩形顶面（高度 120m -> km）
    z_top = z_top_obstacle;
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
for i = 1:size(crowded_zones_km, 1)
    circle = crowded_zones_km(i, :);
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
    
    % 绘制圆形顶面（高度 120m -> km）
    z_top = z_top_obstacle;
    fill3(xCircle, yCircle, z_top * ones(size(theta)), [0.6, 1, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'g', 'LineWidth', 1);
    
    % 绘制圆柱侧面
    [X, Y, Z] = cylinder(radius, 100);
    X = X + cx;
    Y = Y + cy;
    Z = Z * z_top;
    surf(X, Y, Z, 'FaceColor', [0.6, 1, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% ========== 绘制噪音敏感区（矩形，橙色） ==========
for i = 1:size(noise_zones_km, 1)
    rect = noise_zones_km(i, :);
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
    
    % 绘制矩形顶面（高度 120m -> km）
    z_top = z_top_obstacle;
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
depot = dotss_km(21, :);  % 配送中心在第21行（最后一行）
% 配送中心高度设为0（地面）
scatter3(depot(1), depot(2), 0, 200, [0.8, 0.6, 0.8], 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 2);
% 标签位置：在点的上方（y方向），避免遮挡（偏移原为 m，已换 km）
text(depot(1), depot(2) + 200 * m2km, 250 * m2km, '配送中心', 'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% ========== 绘制商家（蓝色圆点，编号1-10） ==========
merchants = dotss_km(1:10, :);  % 商家在第1-10行
% 商家高度设为0（地面）
scatter3(merchants(:, 1), merchants(:, 2), zeros(size(merchants, 1), 1), 100, 'b', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 添加商家标签（标签位置在点的上方y方向，避免遮挡）
for i = 1:size(merchants, 1)
    text(merchants(i, 1), merchants(i, 2) + 200 * m2km, 250 * m2km, ...
        sprintf('商家%d', i), 'FontSize', 8, 'Color', 'b', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% ========== 绘制客户点（绿色圆点，编号11-20） ==========
customers = dotss_km(11:20, :);  % 客户点在第11-20行
% 客户点高度设为0（地面）
scatter3(customers(:, 1), customers(:, 2), zeros(size(customers, 1), 1), 100, 'g', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 添加客户点标签（标签位置在点的上方y方向，避免遮挡）
for i = 1:size(customers, 1)
    text(customers(i, 1), customers(i, 2) + 200 * m2km, 250 * m2km, ...
        sprintf('客户%d', i+10), 'FontSize', 8, 'Color', 'g', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

axisLineWidth = 1.5;
% 设置坐标轴范围（数据已为 km：0--10 km）
xlim([0, 10000 * m2km]);
ylim([0, 10000 * m2km]);
zlim([0, 1000 * m2km]);

grid off;
set(gca, 'Box', 'on', 'LineWidth', axisLineWidth);

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
    
    fprintf('\n========== 各无人机飞行顺序（首尾均为配送中心 21）==========\n');
    % 绘制每架无人机的路径
    for seg_idx = 1:length(segments)
        seg = segments{seg_idx};
        if isempty(seg)
            continue;
        end
        
        % 构建完整路径：配送中心 -> 段内点 -> 配送中心
        path_indices = [21, seg(:)', 21];  % 21是配送中心索引；seg 列向量统一成行再拼
        
        % 控制台输出：该无人机飞行顺序（索引 + 文字）
        fprintf('无人机%d 飞行顺序（节点索引）: ', seg_idx);
        fprintf('%d ', path_indices);
        fprintf('\n');
        lab_parts = cell(1, numel(path_indices));
        for k = 1:numel(path_indices)
            idx = path_indices(k);
            if idx == 21
                lab_parts{k} = '配送中心';
            elseif idx >= 1 && idx <= 10
                lab_parts{k} = sprintf('商家%d', idx);
            elseif idx >= 11 && idx <= 20
                lab_parts{k} = sprintf('客户%d', idx);
            else
                lab_parts{k} = sprintf('节点%d', idx);
            end
        end
        fprintf('  %s\n', strjoin(lab_parts, ' -> '));
        
        % 获取路径点的坐标（与地图一致，使用 km）
        path_coords = zeros(length(path_indices), 3);
        for i = 1:length(path_indices)
            idx = path_indices(i);
            if idx >= 1 && idx <= size(dotss_km, 1)
                path_coords(i, :) = dotss_km(idx, :);
            else
                warning('myapp:indexOutOfRange', '索引 %d 超出范围，跳过', idx);
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
        
        % 添加无人机标签（X 向右偏移，避免与路径重叠）
        if ~isempty(path_coords)
            mid_point = path_coords(ceil(size(path_coords, 1)/2), :);
            text(mid_point(1) + 180 * m2km, mid_point(2), mid_point(3) + 100 * m2km, ...
                sprintf('无人机%d', seg_idx), 'FontSize', 9, ...
                'Color', colors(seg_idx, :), 'FontWeight', 'bold');
        end
    end
    
end

if map_only_env
    % 缩小主图占据区域，右侧与四周留出空白，便于 eastoutside 图例完整落在图窗内
    set(gca, 'Units', 'normalized');
    set(gca, 'Position', [0.10 0.11 0.68 0.80]);
    % 图例：配送中心、商家、客户点、禁飞区、人流密集区、噪音敏感区
    legend_items = {};
    legend_handles = [];
    h1 = scatter3(NaN, NaN, NaN, 200, [0.8, 0.6, 0.8], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    legend_handles = [legend_handles, h1];
    legend_items{end+1} = '配送中心';
    h2 = scatter3(NaN, NaN, NaN, 100, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    legend_handles = [legend_handles, h2];
    legend_items{end+1} = '商家';
    h3 = scatter3(NaN, NaN, NaN, 100, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    legend_handles = [legend_handles, h3];
    legend_items{end+1} = '客户点';
    h4 = fill3(NaN, NaN, NaN, [1, 0, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2);
    legend_handles = [legend_handles, h4];
    legend_items{end+1} = '禁飞区';
    h5 = fill3(NaN, NaN, NaN, [0.6, 1, 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'LineWidth', 1.5);
    legend_handles = [legend_handles, h5];
    legend_items{end+1} = '人流密集区';
    h6 = fill3(NaN, NaN, NaN, [1, 0.5, 0], 'FaceAlpha', 0.4, 'EdgeColor', [1, 0.5, 0], 'LineWidth', 1.5);
    legend_handles = [legend_handles, h6];
    legend_items{end+1} = '噪音敏感区';
    leg = legend(legend_handles, legend_items, 'Location', 'eastoutside', ...
        'Orientation', 'vertical', 'Interpreter', 'none', 'Box', 'on');
    set(leg, 'FontSize', 10);
    % 为底部图例留出边距，导出/保存时不易被裁掉
    drawnow;
    set(figH, 'PaperPositionMode', 'auto');
else
    % 图题：仅显示实验名（如 ablation0416_NSGAIII_G）；Excel 模式且 expName 为空时用数据文件名
    expTitle = strtrim(char(string(expName)));
    if isempty(expTitle)
        [~, expTitle, ~] = fileparts(char(string(excel_file)));
    end
    title(expTitle, 'Interpreter', 'none', 'FontSize', 11);
    set(figH, 'PaperPositionMode', 'auto');
end

hold off;
end

%% ===== Ablation：查找单 run 的 xlsx（多种命名）=====
function p = find_run_excel_file(expDir, expName, runIdx)
% 在 expDir 下尝试多种常见命名；若仍无则通配 *run<idx>*.xlsx
    p = '';
    if nargin < 3 || isempty(runIdx), return; end
    bases = {};
    bases{end+1} = fullfile(expDir, sprintf('%s_run%02d', expName, runIdx));
    bases{end+1} = fullfile(expDir, sprintf('%s_run%d', expName, runIdx));
    bases{end+1} = fullfile(expDir, sprintf('run%02d', runIdx));
    bases{end+1} = fullfile(expDir, sprintf('run%d', runIdx));
    bases{end+1} = fullfile(expDir, sprintf('Run%02d', runIdx));
    bases{end+1} = fullfile(expDir, sprintf('Run%d', runIdx));
    exts = {'.xlsx', '.xls'};
    for b = 1:numel(bases)
        for e = 1:numel(exts)
            cand = [bases{b} exts{e}];
            if isfile(cand)
                p = cand;
                return;
            end
        end
    end
    % 通配兜底（不同机器上命名可能带额外前后缀）
    pat = sprintf('*run%d*.xlsx', runIdx);
    d = dir(fullfile(expDir, pat));
    if ~isempty(d)
        p = fullfile(expDir, d(1).name);
        return;
    end
    pat = sprintf('*Run%d*.xlsx', runIdx);
    d = dir(fullfile(expDir, pat));
    if ~isempty(d)
        p = fullfile(expDir, d(1).name);
    end
end

function [obj_data, dec_data, details_mat] = load_run_from_runs_mat(matPath, runIdx)
% 从 runplatemo_amo_EXP_m 保存的 runs.mat 中取指定 Run 的 Obj_end/Dec_end/Details
    details_mat = [];
    S = load(matPath, 'runs');
    if ~isfield(S, 'runs')
        error('mat 缺少 runs 变量：%s', matPath);
    end
    runs = S.runs;
    runList = arrayfun(@(s)s.Run, runs(:));
    ix = find(runList == runIdx, 1);
    if isempty(ix)
        error('runs.mat 中找不到 Run=%d（现有 Run: %s）', runIdx, mat2str(unique(runList(:))'));
    end
    obj_data = runs(ix).Obj_end;
    dec_data = runs(ix).Dec_end;
    if isfield(runs(ix), 'Details') && ~isempty(runs(ix).Details)
        details_mat = runs(ix).Details;
    end
    if isempty(obj_data) || isempty(dec_data)
        error('Run=%d 的 Obj_end/Dec_end 为空', runIdx);
    end
end

function [runIdx, hv_end] = get_run_index_by_HV_rank_xlsx(summaryPath, hv_rank)
    if exist('readtable', 'file')
        try
            T = readtable(summaryPath, 'Sheet', 'runs', 'VariableNamingRule', 'preserve');
        catch
            T = readtable(summaryPath, 'Sheet', 'runs');
        end
        vars = T.Properties.VariableNames;
        runCol = find(cellfun(@(x) ~isempty(regexpi(x, '^run$')), vars), 1);
        hvCol  = find(cellfun(@(x) ~isempty(regexpi(x, 'HV_end|HVend|HV\\s*end')), vars), 1);
        if isempty(runCol), runCol = 1; end
        if isempty(hvCol),  hvCol  = 2; end
        Run    = T.(vars{runCol});
        HV_end = T.(vars{hvCol});
    else
        raw = xlsread(summaryPath, 'runs');
        if isempty(raw)
            error('runs 表为空：%s', summaryPath);
        end
        Run    = raw(:, 1);
        HV_end = raw(:, 2);
    end
    HV_end = HV_end(:);
    Run    = Run(:);
    valid = isfinite(HV_end) & isfinite(Run);
    HV_end = HV_end(valid);
    Run = Run(valid);

    [~, ord] = sort(HV_end, 'descend');
    idx = min(hv_rank, numel(ord));   % 支持 31 次或更多：按实际条数夹紧
    runIdx = Run(ord(idx));
    hv_end = HV_end(ord(idx));
end

function [runIdx, hv_end] = get_run_index_by_HV_rank_mat(matPath, hv_rank)
    S = load(matPath, 'runs');
    if ~isfield(S, 'runs')
        error('mat 缺少 runs：%s', matPath);
    end
    runs = S.runs;
    if isempty(runs)
        error('runs 为空：%s', matPath);
    end
    HV_end = arrayfun(@(s)s.hv_end, runs(:));
    Run    = arrayfun(@(s)s.Run, runs(:));
    [~, ord] = sort(HV_end, 'descend');
    idx = min(hv_rank, numel(ord));
    runIdx = Run(ord(idx));
    hv_end = HV_end(ord(idx));
end
