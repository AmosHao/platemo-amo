%% ========== 直接运行区域：修改下面的参数后直接运行脚本 ==========
% 取消下面的注释并修改参数即可直接运行画图

% 方式1：使用默认文件（Order_Data/testData/20260214order_obj3minf1.xlsx）
% plot_scatter_from_excel();

% 方式2：指定文件名（会自动在 Order_Data/testData 下查找）
% plot_scatter_from_excel('20260214order_obj3minf1.xlsx');

% 方式3：指定文件名和Sheet
% plot_scatter_from_excel('20260214order_obj3minf1.xlsx', 'obj_data');

% 方式4：完整参数（文件名、Sheet、X列、Y列）
plot_scatter('20260214order.xlsx', 'obj_data', 1, 2);




function plot_scatter(excel_file, sheet_name, x_col, y_col)
% 从指定 Excel 文件读取 Sheet 并绘制散点图
%
% 用法:
%   plot_scatter_from_excel()  
%       使用默认文件 Order_Data/20260214order_obj3minf1.xlsx，Sheet 'obj_data'，第1列vs第2列
%   plot_scatter_from_excel(excel_file)
%       指定 Excel 文件路径，其余默认
%   plot_scatter_from_excel(excel_file, sheet_name)
%       指定文件和 Sheet 名
%   plot_scatter_from_excel(excel_file, sheet_name, x_col, y_col)
%       指定文件、Sheet、X列、Y列
%
% 示例:
%   plot_scatter_from_excel('Order_Data/20260214order_obj3minf1.xlsx');
%   plot_scatter_from_excel('Order_Data/20260214order_obj3minf1.xlsx', 'obj_data', 1, 2);

    % 默认：Excel 文件在 Order_Data/testData 下
    if nargin < 1 || isempty(excel_file)
        script_dir = fileparts(mfilename('fullpath'));
        excel_file = fullfile(script_dir, '..', 'Order_Data', 'testData', '20260214order_obj3minf1.xlsx');
    else
        % 若只给出文件名（无路径），则在 Order_Data/testData 下查找
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
    
    if ~isfile(excel_file)
        error('文件不存在: %s', excel_file);
    end
    
    % 读取指定 Sheet
    data = readmatrix(excel_file, 'Sheet', sheet_name);
    
    if isempty(data)
        error('Sheet "%s" 为空或读取失败', sheet_name);
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
    figure('Name', 'Scatter from Excel');
    scatter(x, y, 36, 'b', 'filled');
    xlabel(sprintf('第%d列', x_col));
    ylabel(sprintf('第%d列', y_col));
    title(sprintf('散点图: %s (Sheet: %s)', get_filename_only(excel_file), sheet_name));
    grid on;
    
    % 若为 obj_data 且为前两列，使用常用标签
    if strcmpi(sheet_name, 'obj_data') && x_col == 1 && y_col == 2
        xlabel('f1 (总能量/J)');
        ylabel('f2 (总时间/s)');
        title(sprintf('Pareto 前沿: %s', get_filename_only(excel_file)));
    end
end

function name = get_filename_only(path)
    [~, name, ext] = fileparts(path);
    name = [name ext];
end


