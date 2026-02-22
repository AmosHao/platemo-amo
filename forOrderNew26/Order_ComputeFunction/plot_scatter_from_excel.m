%% ========== 直接运行区域：修改下面的参数后直接运行脚本 ==========
% 取消下面的注释并修改参数即可直接运行画图

% 方式1：使用默认文件（Order_Data/testData/20260214order_obj3minf1.xlsx）
% plot_scatter_from_excel();

% 方式2：指定文件名（会自动在 Order_Data/testData 下查找）
% plot_scatter_from_excel('20260214order_obj3minf1.xlsx');

% 方式3：指定文件名和Sheet
% plot_scatter_from_excel('20260214order_obj3minf1.xlsx', 'obj_data');

% 方式4：完整参数（文件名、Sheet、X列、Y列）
% plot_scatter('20260222order.xlsx', 'obj_data', 1, 2);

% 方式5：完整参数 + 筛选可行解（filter_feasible=true 时只保留 details_data 中 Feasible=1 的解）
plot_scatter('2026022221order.xlsx', 'obj_data', 1, 2, true);




function plot_scatter(excel_file, sheet_name, x_col, y_col, filter_feasible)
% 从指定 Excel 文件读取 Sheet 并绘制散点图
%
% 用法:
%   plot_scatter(excel_file)
%       指定 Excel 文件路径，其余默认
%   plot_scatter(excel_file, sheet_name)
%       指定文件和 Sheet 名
%   plot_scatter(excel_file, sheet_name, x_col, y_col)
%       指定文件、Sheet、X列、Y列
%   plot_scatter(excel_file, sheet_name, x_col, y_col, filter_feasible)
%       filter_feasible=true 时，先读取 details_data Sheet，筛出 Feasible=1 的行索引，
%       再从 sheet_name 中只取这些行绘图（需要 Excel 中存在 details_data Sheet）
%
% 示例:
%   plot_scatter('20260222order.xlsx', 'obj_data', 1, 2, true);   % 仅绘制可行解

    % 默认：Excel 文件在 Order_Data/testData 下
    if nargin < 1 || isempty(excel_file)
        script_dir = fileparts(mfilename('fullpath'));
        excel_file = fullfile(script_dir, '..', 'Order_Data', 'testData', '20260214order_obj3minf1.xlsx');
    else
        [filepath, ~, ~] = fileparts(excel_file);
        if isempty(filepath)
            script_dir = fileparts(mfilename('fullpath'));
            excel_file = fullfile(script_dir, '..', 'Order_Data', 'testData', excel_file);
        end
    end
    
    if nargin < 2 || isempty(sheet_name)
        sheet_name = 'obj_data';
    end
    if nargin < 3 || isempty(x_col)
        x_col = 1;
    end
    if nargin < 4 || isempty(y_col)
        y_col = 2;
    end
    if nargin < 5 || isempty(filter_feasible)
        filter_feasible = false;
    end
    
    if ~isfile(excel_file)
        error('文件不存在: %s', excel_file);
    end
    
    % ---- 筛选可行解 ----
    % 读取 details_data，找出 Feasible=1（最后一列）的行索引
    feasible_rows = [];
    if filter_feasible
        try
            details = readmatrix(excel_file, 'Sheet', 'details_data');
            if ~isempty(details)
                feasible_col = details(:, end);   % 最后一列为 Feasible 标志
                feasible_rows = find(feasible_col == 1);
                fprintf('[筛选] details_data 共 %d 行，其中可行解 %d 个\n', ...
                    size(details, 1), numel(feasible_rows));
                if isempty(feasible_rows)
                    warning('myapp:noFeasible', '%s', '没有找到任何可行解（Feasible=1），将显示所有解');
                    filter_feasible = false;
                end
            else
                warning('myapp:emptyDetails', '%s', 'details_data Sheet 为空，忽略可行解筛选');
                filter_feasible = false;
            end
        catch e
            warning('myapp:readDetailsFailed', '读取 details_data 失败（%s），忽略可行解筛选', e.message);
            filter_feasible = false;
        end
    end
    
    % 读取目标 Sheet
    data = readmatrix(excel_file, 'Sheet', sheet_name);
    
    if isempty(data)
        error('Sheet "%s" 为空或读取失败', sheet_name);
    end
    
    % 若已筛出可行行索引，先过滤数据行
    if filter_feasible && ~isempty(feasible_rows)
        valid_in_range = feasible_rows(feasible_rows <= size(data, 1));
        data = data(valid_in_range, :);
    end
    
    % 取 X、Y 列（去掉可能的 NaN 行）
    x = data(:, x_col);
    y = data(:, y_col);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    
    if isempty(x)
        error('有效数据为空，请检查列索引与数据');
    end
    
    % 绘制散点图
    fig_name = 'Scatter from Excel';
    if filter_feasible
        fig_name = 'Scatter from Excel (仅可行解)';
    end
    figure('Name', fig_name);
    scatter(x, y, 36, 'b', 'filled');
    xlabel(sprintf('第%d列', x_col));
    ylabel(sprintf('第%d列', y_col));
    title(sprintf('散点图: %s (Sheet: %s)', get_filename_only(excel_file), sheet_name));
    grid on;
    
    % 若为 obj_data 且为前两列，使用常用标签
    if strcmpi(sheet_name, 'obj_data') && x_col == 1 && y_col == 2
        xlabel('f1 (总能量/J)');
        ylabel('f2 (总时间/s)');
        feasible_tag = '';
        if filter_feasible
            feasible_tag = sprintf(' [仅可行解 %d 个]', numel(x));
        end
        title(sprintf('Pareto 前沿: %s%s', get_filename_only(excel_file), feasible_tag));
    end
end

function name = get_filename_only(path)
    [~, name, ext] = fileparts(path);
    name = [name ext];
end


