% =====================================================================
% 按 30 次运行的最终 HV 排名，画出指定名次那次运行的目标值 (f1-f2) 散点图
% 数据：Order_Data/ablation/<实验名>/ 下 summary 的 runs 表（HV_end 排名）+ 对应 runXX.xlsx 的 obj_data
% =====================================================================
% 用法：在编辑器中按 F5 运行；或命令窗口：
%   plot_ablation_obj_by_HV_rank(expNames, hv_rank)
% 参数：
%   expNames - 实验名（与 ablation 下子文件夹一致），cell 或字符。留空则自动扫描所有实验
%   hv_rank  - 要画的排名，1=HV 最高的一次，2=第 2 名，… 30=第 30 名。默认 2
% 输出：Order_Data/ablation/all_HVrank<k>_obj_scatter.png（及 .fig），所有实验在一张图
% =====================================================================

% ---------- 直接运行本文件时执行（F5）：可修改下面两行 ----------
expNames = {'ablation06_A_amos', 'ablation06_A_amos_noCluster','ablation06_A_amos_allBYJC'};   % 留空=扫描所有实验；或指定如 {'ablation_A_amos'}
hv_rank  = 15;    % 30 次中按最终 HV 排名，画第几名的那次（1=最好，2=第2名…）
plot_ablation_obj_by_HV_rank_inside(expNames, hv_rank);

function plot_ablation_obj_by_HV_rank_inside(expNames, hv_rank)
    % 与文件名一致的对入口，仅转发到内部实现（内部函数名不与文件名一致）
    run_plot_ablation_obj_by_HV_rank(expNames, hv_rank);
end

function run_plot_ablation_obj_by_HV_rank(expNames, hv_rank)
    if nargin < 1
        expNames = {};
    end
    if nargin < 2
        hv_rank = 2;
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
    end

    hv_rank = max(1, min(30, round(hv_rank)));

    % 每个实验一种颜色，画在一张图上
    colors = lines(max(length(expNames), 1));

    figure('Name', sprintf('各实验 HV排名第%d 目标值', hv_rank), 'NumberTitle', 'off');
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
            runIdx = get_run_index_by_HV_rank(summaryPath, hv_rank);
            if isempty(runIdx)
                warning('实验 %s 无法得到 HV 排名第 %d 的 run 序号，跳过。', expName, hv_rank);
                continue;
            end
            runPath = fullfile(ablation_dir, expName, [expName '_run' sprintf('%02d', runIdx) '.xlsx']);
            if ~isfile(runPath)
                warning('未找到 run 文件: %s，跳过。', runPath);
                continue;
            end
            Obj = read_obj_data(runPath);
            runPath = fullfile(ablation_dir, expName, [expName '_run' sprintf('%02d', runIdx) '.xlsx']);
            if ~isfile(runPath)
                warning('未找到 run 文件: %s，跳过。', runPath);
                continue;
            end
            Obj = read_obj_data(runPath);
            if isempty(Obj)
                warning('实验 %s run%02d 的 obj_data 为空，跳过。', expName, runIdx);
                continue;
            end
            f1 = Obj(:, 1);
            f2 = Obj(:, 2);
            scatter(f1, f2, 24, colors(i,:), 'filled', 'DisplayName', strrep(expName, '_', '\_'));
            plottedNames{end+1} = expName; %#ok<AGROW>
        catch me
            warning('处理实验 %s 时出错: %s', expName, me.message);
        end
    end
    hold off;
    xlabel('f1 (能耗)', 'FontSize', 11);
    ylabel('f2 (时间)', 'FontSize', 11);
    title(sprintf('各实验 30 次中 HV 排名第 %d 的那次 — 目标值', hv_rank));
    if ~isempty(plottedNames)
        legend(plottedNames, 'Interpreter', 'none', 'Location', 'best');
    end
    grid on;
    set(gca, 'FontSize', 10);

    if ~isempty(plottedNames)
        figName = fullfile(ablation_dir, ['all_HVrank' num2str(hv_rank) '_obj_scatter']);
        saveas(gcf, [figName '.png']);
        savefig(gcf, [figName '.fig']);
        fprintf('已保存: %s.png / .fig（HV 排名第 %d，共 %d 个实验）\n', figName, hv_rank, length(plottedNames));
    end
end

function runIdx = get_run_index_by_HV_rank(summaryPath, hv_rank)
    % 读 summary 的 runs 表：Run, HV_end, IGD_end；按 HV_end 降序排名，取第 hv_rank 名对应的 Run
    if exist('readtable', 'file')
        try
            T = readtable(summaryPath, 'Sheet', 'runs', 'VariableNamingRule', 'preserve');
        catch
            T = readtable(summaryPath, 'Sheet', 'runs');
        end
        vars = T.Properties.VariableNames;
        runCol = find(cellfun(@(x) ~isempty(regexpi(x, '^run$')), vars), 1);
        hvCol  = find(cellfun(@(x) ~isempty(regexpi(x, 'HV_end|HVend')), vars), 1);
        if isempty(runCol), runCol = 1; end
        if isempty(hvCol),  hvCol  = 2; end
        Run    = T.(vars{runCol});
        HV_end = T.(vars{hvCol});
    else
        raw = xlsread(summaryPath, 'runs');
        if isempty(raw)
            runIdx = [];
            return;
        end
        Run    = raw(:, 1);
        HV_end = raw(:, 2);
    end
    [~, ord] = sort(HV_end, 'descend');
    idx = min(hv_rank, numel(ord));
    runIdx = Run(ord(idx));
end

function Obj = read_obj_data(runPath)
    if exist('readmatrix', 'file')
        Obj = readmatrix(runPath, 'Sheet', 'obj_data');
    else
        Obj = xlsread(runPath, 'obj_data');
    end
    if isempty(Obj) || size(Obj, 2) < 2
        Obj = [];
    end
end
