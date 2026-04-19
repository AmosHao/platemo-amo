% =====================================================================
% 根据 runplatemo_amo_EXP.m 生成的消融实验总数据表，绘制 30 次运行平均 HV 折线图
% 数据来源：Order_Data/ablation/<实验名>/<实验名>_summary.xlsx 中的 hv_mean_curve 表
% =====================================================================
% 用法：在编辑器中直接按 F5 运行本文件即可（下方可改 expNames）。
%   也可在命令窗口调用：plot_ablation_HV_mean() 或 plot_ablation_HV_mean({'ablation_A_amos'})
% 输出：图片保存到 Order_Data/ablation/all_HV_mean.png（及 .fig），所有实验在一张图
% =====================================================================

% ---------- 直接运行本文件时执行（F5）：写死一个绝对路径即可 ----------
% 这个目录下应包含子文件夹：ablation0415_NSGAIII、ablation0415_NSGAIII_G、...
ablation_dir = 'D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\newdata2604\ablation0416';
% algorithms   = {'NSGAIII_ori','NSGAIII_G','NSGAIII_MG','NSGAIII_MZG'};  % 只画这些算法；留空则扫描该目录下全部实验
% algorithms   = {'NSGAIII_ori','NSGAIII_G','noJF_NSGAIII_G'}; 
algorithms = {
    % 'NSGAIII_ori','NSGAIII_G',...
    % 'PeEA_ori','PeEA_G',...
    % 'NSGAII_ori','NSGAII_G',...
    % 'CLIA_ori','CLIA_G',...
    % 'SMSEMOA_ori','SMSEMOA_G',...
    % 'MOEAD_ori','MOEAD_G',...
    'MOEADD_ori','MOEADD_G'
};
expNames     = {};  % 一般不用填：会由 algorithms 自动生成；或手填具体实验文件夹名
plot_ablation_HV_mean_inside(expNames, ablation_dir, algorithms);

function plot_ablation_HV_mean_inside(expNames, ablation_dir, algorithms)
    if nargin < 1 || isempty(expNames), expNames = {}; end
    if nargin < 2, ablation_dir = ''; end
    if nargin < 3, algorithms = {}; end
    if ischar(expNames) || isstring(expNames)
        expNames = cellstr(expNames);
    end
    if ischar(algorithms) || isstring(algorithms)
        algorithms = cellstr(algorithms);
    end

    % 只要传了 ablation_dir（建议绝对路径），就只在该目录下读；不存在则直接报错
    if isempty(ablation_dir)
        script_dir   = fileparts(mfilename('fullpath'));
        ablation_dir = fullfile(script_dir, '..', 'Order_Data', 'ablation');
    end
    if ~isfolder(ablation_dir)
        error('Order_EXP_Function:NoAblationDir', '未找到目录: %s', ablation_dir);
    end
    expNames = resolve_expNames_from_algorithms_or_keep(ablation_dir, expNames, algorithms);

    if isempty(expNames)
        if ~isfolder(ablation_dir)
            error('Order_EXP_Function:NoAblationDir', '未找到目录: %s', ablation_dir);
        end
        list = dir(ablation_dir);
        expNames = {};
        for k = 1:length(list)
            if ~list(k).isdir || strcmp(list(k).name, '.') || strcmp(list(k).name, '..')
                continue;
            end
            folderName = list(k).name;
            summaryPath = fullfile(ablation_dir, folderName, [folderName '_summary.xlsx']);
            matPath = fullfile(ablation_dir, folderName, [folderName '_runs.mat']);
            if isfile(summaryPath) || isfile(matPath)
                expNames{end+1} = folderName; %#ok<AGROW>
            end
        end
        if isempty(expNames)
            fprintf('目录下未找到任何带 _summary.xlsx 或 _runs.mat 的子文件夹：%s\n', ablation_dir);
            return;
        end
        fprintf('自动发现 %d 个实验: %s\n', length(expNames), strjoin(expNames, ', '));
    end

    % 每个实验一种颜色；第 1 条空心方块、第 2 条空心圆，其余按 markers 循环
    % prepend_corner：true 时在首点前插入点 (0, yLo)，yLo 为即将设置的 y 轴下界（图左下角），不用 (0,0)；false 仅画表内数据
    prepend_corner = true;
    lineColors = {'b', 'r', 'g', 'm', 'c', 'k', [0.85 0.33 0.1], [0.47 0.67 0.19]};
    markers    = {'^', 'd', 'v', '*', 'p', 'h', '^', 'd'};  % 第 3 条起
    markerSize = 8;
    fontName   = 'Times New Roman';
    axisLineWidth = 1.5;  % 坐标轴线宽（四边）

    nExp = numel(expNames);
    allGen = cell(nExp, 1);
    allHV  = cell(nExp, 1);
    rowOk  = false(nExp, 1);
    for ii = 1:nExp
        try
            [g, h] = read_hv_mean_curve_any(ablation_dir, expNames{ii});
            allGen{ii} = g(:);
            allHV{ii}  = h(:);
            rowOk(ii)  = true;
        catch me
            warning('读取实验 %s 时出错: %s', expNames{ii}, me.message);
        end
    end
    if ~any(rowOk)
        return;
    end

    hvMinData = inf;
    hvMaxData = -inf;
    genMaxData = 0;
    for ii = find(rowOk).'
        g = allGen{ii};
        h = allHV{ii};
        v = isfinite(g) & isfinite(h);
        if ~any(v)
            continue;
        end
        hvMinData = min(hvMinData, min(h(v)));
        hvMaxData = max(hvMaxData, max(h(v)));
        genMaxData = max(genMaxData, max(g(v)));
    end

    if hvMaxData > hvMinData
        span = hvMaxData - hvMinData;
        padY = max(0.02 * span, 1e-9);
        yLo = hvMinData - padY;
        yHi = hvMaxData + padY;
    else
        padY = max(0.001 * max(abs(hvMaxData), 1), 1e-9);
        yLo = hvMaxData - padY;
        yHi = hvMaxData + padY;
    end

    figure('Name', '各实验 30 次平均 HV', 'NumberTitle', 'off');
    hold on;
    plottedNames = {};
    for i = 1:length(expNames)
        if ~rowOk(i)
            continue;
        end
        expName = expNames{i};
        Gen     = allGen{i};
        HV_mean = allHV{i};
        if prepend_corner
            Gen     = [0; Gen(:)];
            HV_mean = [yLo; HV_mean(:)];
        end

        c = lineColors{1 + mod(i - 1, length(lineColors))};
        if i == 1
            mk = 's';
        elseif i == 2
            mk = 'o';
        else
            mk = markers{1 + mod(i - 3, length(markers))};
        end
        if ischar(c) || isstring(c)
            lineColor = char(c);
            lineColor = lineColor(1);
        else
            lineColor = c;
        end
        plot(Gen, HV_mean, '-', 'Color', lineColor, 'Marker', mk, ...
            'LineWidth', 1.5, 'MarkerSize', markerSize, ...
            'MarkerFaceColor', 'none', 'MarkerEdgeColor', lineColor, ...
            'DisplayName', strrep(expName, '_', '\_'));
        plottedNames{end+1} = expName; %#ok<AGROW>
    end
    hold off;
    xlabel('Number of Function Evaluations', 'FontSize', 11, 'FontName', fontName);
    ylabel('HV', 'FontSize', 11, 'FontName', fontName);
    leg = legend(plottedNames, 'Interpreter', 'none', 'Location', 'best');
    set(leg, 'FontName', fontName, 'FontSize', 10);
    grid off;
    set(gca, 'FontSize', 10, 'FontName', fontName, 'Box', 'on', 'LineWidth', axisLineWidth);
    if prepend_corner && isfinite(genMaxData) && genMaxData >= 0
        xlim([0, genMaxData * 1.02 + eps]);
    end
    if isfinite(yLo) && isfinite(yHi)
        ylim([yLo, yHi]);
    end

    if isempty(plottedNames)
        return;
    end
    figName = fullfile(ablation_dir, 'all_HV_mean');
    saveas(gcf, [figName '.png']);
    savefig(gcf, [figName '.fig']);
    fprintf('已保存: %s.png / .fig\n', figName);
end

function [Gen, HV_mean] = read_hv_mean_curve_xlsx(summaryPath)
    % 读取 hv_mean_curve：表头 Gen, HV_mean, HV_std，数据从 A2 开始
    if exist('readtable', 'file')
        try
            T = readtable(summaryPath, 'Sheet', 'hv_mean_curve', 'VariableNamingRule', 'preserve');
        catch
            T = readtable(summaryPath, 'Sheet', 'hv_mean_curve');
        end
        vars = T.Properties.VariableNames;
        genCol  = find(cellfun(@(x) ~isempty(regexpi(x, '^gen$')), vars), 1);
        meanCol = find(cellfun(@(x) ~isempty(regexpi(x, 'HV_mean|HVmean')), vars), 1);
        if isempty(genCol),  genCol = 1; end
        if isempty(meanCol), meanCol = 2; end
        Gen     = T.(vars{genCol});
        HV_mean = T.(vars{meanCol});
    else
        raw = xlsread(summaryPath, 'hv_mean_curve');
        if isempty(raw)
            error('hv_mean_curve 表为空或读取失败。');
        end
        Gen     = raw(:, 1);
        HV_mean = raw(:, 2);
    end
    valid = isfinite(Gen) & isfinite(HV_mean);
    Gen     = Gen(valid);
    HV_mean = HV_mean(valid);
end

function [Gen, HV_mean] = read_hv_mean_curve_mat(matPath)
    S = load(matPath, 'summary');
    if ~isfield(S, 'summary') || ~isfield(S.summary, 'hv_mean')
        error('mat 文件缺少 summary.hv_mean：%s', matPath);
    end
    HV_mean = S.summary.hv_mean(:);
    Gen = (1:numel(HV_mean))';
    valid = isfinite(Gen) & isfinite(HV_mean);
    Gen = Gen(valid);
    HV_mean = HV_mean(valid);
end

function [Gen, HV_mean] = read_hv_mean_curve_any(ablation_dir, expName)
    summaryPath = fullfile(ablation_dir, expName, [expName '_summary.xlsx']);
    matPath     = fullfile(ablation_dir, expName, [expName '_runs.mat']);
    if isfile(summaryPath)
        [Gen, HV_mean] = read_hv_mean_curve_xlsx(summaryPath);
    elseif isfile(matPath)
        [Gen, HV_mean] = read_hv_mean_curve_mat(matPath);
    else
        error('未找到总数据（xlsx 或 mat）：%s / %s', summaryPath, matPath);
    end
end

function expNames = resolve_expNames_from_algorithms_or_keep(ablation_dir, expNames, algorithms)
    if ~isempty(algorithms)
        [~, prefix] = fileparts(ablation_dir);
        expNames = cell(size(algorithms));
        for i = 1:numel(algorithms)
            expNames{i} = [prefix '_' char(algorithms{i})];
        end
        return;
    end
    if ischar(expNames) || isstring(expNames)
        expNames = cellstr(expNames);
    end
end
