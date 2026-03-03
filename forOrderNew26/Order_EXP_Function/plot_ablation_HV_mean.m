% =====================================================================
% 根据 runplatemo_amo_EXP.m 生成的消融实验总数据表，绘制 30 次运行平均 HV 折线图
% 数据来源：Order_Data/ablation/<实验名>/<实验名>_summary.xlsx 中的 hv_mean_curve 表
% =====================================================================
% 用法：在编辑器中直接按 F5 运行本文件即可（下方可改 expNames）。
%   也可在命令窗口调用：plot_ablation_HV_mean() 或 plot_ablation_HV_mean({'ablation_A_amos'})
% 输出：图片保存到 Order_Data/ablation/all_HV_mean.png（及 .fig），所有实验在一张图
% =====================================================================

% ---------- 直接运行本文件时执行（F5）：可修改 expNames 指定要画图的实验 ----------
expNames = {'ablation06_A_amos', 'ablation06_A_amos_noCluster','ablation06_A_amos_allBYJC'};   % 留空 = 自动扫描 ablation 下所有实验；或指定如 {'ablation_A_amos', 'ablation_A_amos_noCluster'}
plot_ablation_HV_mean_inside(expNames);

function plot_ablation_HV_mean_inside(expNames)
    if nargin < 1
        expNames = {};
    end
    if ischar(expNames) || isstring(expNames)
        expNames = cellstr(expNames);
    end

    script_dir   = fileparts(mfilename('fullpath'));
    ablation_dir = fullfile(script_dir, '..', 'Order_Data', 'ablation');

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
            if isfile(summaryPath)
                expNames{end+1} = folderName; %#ok<AGROW>
            end
        end
        if isempty(expNames)
            fprintf('ablation 下未找到任何带 _summary.xlsx 的子文件夹。\n');
            return;
        end
        fprintf('自动发现 %d 个实验: %s\n', length(expNames), strjoin(expNames, ', '));
    end

    % 每个实验一种线型（颜色+标记）
    styles = {'b-o', 'r-s', 'g-^', 'm-d', 'c-v', 'k-*', [0.85 0.33 0.1], [0.47 0.67 0.19]};
    markerSize = 4;

    figure('Name', '各实验 30 次平均 HV', 'NumberTitle', 'off');
    hold on;
    plottedNames = {};
    for i = 1:length(expNames)
        expName = expNames{i};
        summaryPath = fullfile(ablation_dir, expName, [expName '_summary.xlsx']);
        if ~isfile(summaryPath)
            warning('未找到总数据表: %s，跳过。', summaryPath);
            continue;
        end
        try
            [Gen, HV_mean] = read_hv_mean_curve(summaryPath);
        catch me
            warning('读取实验 %s 时出错: %s', expName, me.message);
            continue;
        end
        sty = styles{1 + mod(i - 1, length(styles))};
        if ischar(sty) || isstring(sty)
            plot(Gen, HV_mean, sty, 'LineWidth', 1.5, 'MarkerSize', markerSize, ...
                'MarkerFaceColor', sty(1), 'DisplayName', strrep(expName, '_', '\_'));
        else
            plot(Gen, HV_mean, '-o', 'Color', sty, 'LineWidth', 1.5, 'MarkerSize', markerSize, ...
                'MarkerFaceColor', sty, 'DisplayName', strrep(expName, '_', '\_'));
        end
        plottedNames{end+1} = expName; %#ok<AGROW>
    end
    hold off;
    xlabel('代数 (Gen)', 'FontSize', 11);
    ylabel('HV (30次平均)', 'FontSize', 11);
    title('各实验 30 次运行平均 HV');
    legend(plottedNames, 'Interpreter', 'none', 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 10);

    if isempty(plottedNames)
        return;
    end
    figName = fullfile(ablation_dir, 'all_HV_mean');
    saveas(gcf, [figName '.png']);
    savefig(gcf, [figName '.fig']);
    fprintf('已保存: %s.png / .fig\n', figName);
end

function [Gen, HV_mean] = read_hv_mean_curve(summaryPath)
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
