%% ========== 直接运行区域 ==========
% 从 Excel 的 clust_data Sheet 读取数据并绘制聚类着色散点图
% 取消下方注释后直接运行：
%
% plot_clust_scatter_from_excel('2026022224_03order.xlsx');
% plot_clust_scatter_from_excel('2026022224_03order.xlsx', true);  % 仅可行解

function plot_clust_scatter(obj_matrix, clust_tags, title_str)
% 按聚类标签对 Pareto 前沿解着色绘制散点图
%
% 用法（直接传矩阵）：
%   plot_clust_scatter(obj_matrix, clust_tags)
%   plot_clust_scatter(obj_matrix, clust_tags, title_str)
%
% 输入:
%   obj_matrix  - N×2 目标值矩阵 [f1, f2]
%   clust_tags  - N×1 聚类编号向量（整数，每解对应一个类号）
%   title_str   - （可选）图标题字符串

    if nargin < 3 || isempty(title_str)
        title_str = 'Pareto 前沿聚类分布';
    end

    % 去掉无效行（inf / nan）
    valid = isfinite(obj_matrix(:,1)) & isfinite(obj_matrix(:,2)) & isfinite(clust_tags);
    obj_matrix  = obj_matrix(valid, :);
    clust_tags  = clust_tags(valid);

    if isempty(obj_matrix)
        warning('plot_clust_scatter: 没有有效数据可绘制');
        return;
    end

    % 对聚类编号重新紧凑编号（1,2,3,...,K）
    unique_tags = unique(clust_tags);
    K = numel(unique_tags);
    compact_tags = zeros(size(clust_tags));
    for k = 1:K
        compact_tags(clust_tags == unique_tags(k)) = k;
    end

    % 生成 K 种区分度高的颜色（使用 HSV 色系均匀分布）
    if K <= 12
        cmap = lines(K);  % MATLAB 内置高对比色表
    else
        cmap = hsv(K);
    end

    figure('Name', ['聚类散点图: ' title_str], 'Color', 'w');
    hold on;

    legend_entries = cell(K, 1);
    for k = 1:K
        idx = compact_tags == k;
        scatter(obj_matrix(idx, 1), obj_matrix(idx, 2), 48, ...
            cmap(k, :), 'filled', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.85);
        legend_entries{k} = sprintf('Cluster %d (%d 解)', unique_tags(k), sum(idx));
    end

    % 标注每类的聚类中心
    for k = 1:K
        idx = compact_tags == k;
        cx = mean(obj_matrix(idx, 1));
        cy = mean(obj_matrix(idx, 2));
        plot(cx, cy, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
    end

    xlabel('f_1 (总能量 / J)', 'FontSize', 12);
    ylabel('f_2 (总时间 / s)', 'FontSize', 12);
    title(sprintf('Pareto 前沿聚类分布: %s\n（共 %d 类，%d 个解）', ...
        title_str, K, size(obj_matrix, 1)), 'FontSize', 12);
    legend(legend_entries, 'Location', 'best', 'FontSize', 9);
    grid on;
    box on;
    hold off;

    fprintf('[聚类散点图] 共 %d 类，%d 个解\n', K, size(obj_matrix, 1));
    for k = 1:K
        idx = compact_tags == k;
        fprintf('  Cluster %2d: %3d 解  f1=[%.1f, %.1f]  f2=[%.1f, %.1f]\n', ...
            unique_tags(k), sum(idx), ...
            min(obj_matrix(idx,1)), max(obj_matrix(idx,1)), ...
            min(obj_matrix(idx,2)), max(obj_matrix(idx,2)));
    end
end


function plot_clust_scatter_from_excel(excel_file, filter_feasible)
% 从 Excel 文件的 clust_data Sheet 读取数据，绘制聚类着色散点图
%
% 用法:
%   plot_clust_scatter_from_excel('2026022224_03order.xlsx')
%   plot_clust_scatter_from_excel('2026022224_03order.xlsx', true)  % 仅可行解

    if nargin < 1 || isempty(excel_file)
        script_dir = fileparts(mfilename('fullpath'));
        excel_file = fullfile(script_dir, '..', 'Order_Data', 'testData', '2026022224_03order.xlsx');
    else
        [filepath, ~, ~] = fileparts(excel_file);
        if isempty(filepath)
            script_dir = fileparts(mfilename('fullpath'));
            excel_file = fullfile(script_dir, '..', 'Order_Data', 'testData', excel_file);
        end
    end
    if nargin < 2 || isempty(filter_feasible)
        filter_feasible = false;
    end

    if ~isfile(excel_file)
        error('文件不存在: %s', excel_file);
    end

    % 读取 clust_data（第1行为表头）
    raw = readcell(excel_file, 'Sheet', 'clust_data');
    if size(raw, 1) < 2
        error('clust_data Sheet 数据不足（需至少1行表头 + 1行数据）');
    end
    data = cell2mat(raw(2:end, :));  % 跳过表头行

    % 可行解筛选
    if filter_feasible
        try
            details = readmatrix(excel_file, 'Sheet', 'details_data');
            feasible_rows = find(details(:, end) == 1);
            if isempty(feasible_rows)
                warning('myapp:noFeasible', '没有找到可行解，将显示所有解');
            else
                valid_in_range = feasible_rows(feasible_rows <= size(data, 1));
                data = data(valid_in_range, :);
                fprintf('[筛选] 仅显示 %d 个可行解\n', size(data, 1));
            end
        catch e
            warning('myapp:readFailed', '读取 details_data 失败: %s', e.message);
        end
    end

    obj_matrix  = data(:, 1:2);
    clust_tags  = data(:, 3);

    [~, fname, ext] = fileparts(excel_file);
    plot_clust_scatter(obj_matrix, clust_tags, [fname ext]);
end
