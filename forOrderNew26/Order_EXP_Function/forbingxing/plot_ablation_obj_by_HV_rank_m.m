% =====================================================================
% 按 30 次运行的最终 HV 排名，画出指定名次那次运行的目标值 (f1-f2) 散点图
% 数据：优先读 Order_Data/ablation/<实验名>/ 下的 _summary.xlsx（旧格式）
%      若不存在，则读 <实验名>_runs.mat（并行 mat 格式）
% =====================================================================
% 用法：在编辑器中按 F5 运行；或命令窗口：
%   plot_ablation_obj_by_HV_rank(expNames, hv_rank)
% 参数：
%   expNames - 实验名（与 ablation 下子文件夹一致），cell 或字符。留空则自动扫描所有实验
%   hv_rank  - 要画的排名，1=HV 最高的一次，2=第 2 名，… 30=第 30 名。默认 2
% 输出：Order_Data/ablation/all_HVrank<k>_obj_scatter.png（及 .fig），所有实验在一张图
% =====================================================================

% ---------- 直接运行本文件时执行（F5）：写死一个绝对路径即可 ----------
% 这个目录下应包含子文件夹：ablation0415_NSGAIII、ablation0415_NSGAIII_G、...
ablation_dir = 'D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\newdata2604\ablation0416';
% algorithms   = {'NSGAIII_ori','NSGAIII_G','NSGAIII_MG','NSGAIII_MZG'};  % 只画这些算法；留空则扫描该目录下全部实验
% algorithms   = {'NSGAIII_ori','NSGAIII_G'}; 
algorithms = {
    'NSGAIII_ori','NSGAIII_G',...
    % 'PeEA_ori','PeEA_G',...
    % 'NSGAII_ori','NSGAII_G',...
    % 'CLIA_ori','CLIA_G',...
    % 'SMSEMOA_ori','SMSEMOA_G',...
    % 'MOEAD_ori','MOEAD_G',...
    % 'MOEADD_ori','MOEADD_G'
};
hv_rank      = 16;
expNames     = {};   % 一般不用填：会由 algorithms 自动生成；或手填具体实验文件夹名
plot_ablation_obj_by_HV_rank_inside(expNames, hv_rank, ablation_dir, algorithms);

function plot_ablation_obj_by_HV_rank_inside(expNames, hv_rank, ablation_dir, algorithms)
    % 与文件名一致的对入口，仅转发到内部实现（内部函数名不与文件名一致）
    run_plot_ablation_obj_by_HV_rank(expNames, hv_rank, ablation_dir, algorithms);
end

function run_plot_ablation_obj_by_HV_rank(expNames, hv_rank, ablation_dir, algorithms)
    if nargin < 1
        expNames = {};
    end
    if nargin < 2
        hv_rank = 2;
    end
    if nargin < 3
        ablation_dir = '';
    end
    if nargin < 4
        algorithms = {};
    end
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
    end

    hv_rank = max(1, min(30, round(hv_rank)));

    fontName        = 'Times New Roman';
    axisLineWidth   = 1.5;

    % 每个实验一种颜色，画在一张图上
    colors = lines(max(length(expNames), 1));

    figure('Name', sprintf('各实验 HV排名第%d 目标值', hv_rank), 'NumberTitle', 'off');
    hold on;
    plottedNames = {};
    for i = 1:length(expNames)
        expName = expNames{i};
        summaryPath = fullfile(ablation_dir, expName, [expName '_summary.xlsx']);
        matPath     = fullfile(ablation_dir, expName, [expName '_runs.mat']);
        hasXlsx     = isfile(summaryPath);
        hasMat      = isfile(matPath);
        if ~hasXlsx && ~hasMat
            warning('未找到总数据（xlsx 或 mat）: %s / %s，跳过。', summaryPath, matPath);
            continue;
        end
        try
            if hasXlsx
                runIdx = get_run_index_by_HV_rank_xlsx(summaryPath, hv_rank);
                Obj = read_obj_data_xlsx(fullfile(ablation_dir, expName, [expName '_run' sprintf('%02d', runIdx) '.xlsx']));
            else
                [runIdx, Obj] = get_obj_by_HV_rank_mat(matPath, hv_rank);
            end
            if isempty(runIdx)
                warning('实验 %s 无法得到 HV 排名第 %d 的 run 序号，跳过。', expName, hv_rank);
                continue;
            end
            if isempty(Obj)
                warning('实验 %s run%02d 的 obj_data 为空，跳过。', expName, runIdx);
                continue;
            end
            f1 = Obj(:, 1);
            f2 = Obj(:, 2);
            ptArea = 36;
            if i == 1
                mk  = 's';
                mfc = 'b';
            elseif i == 2
                mk  = 'o';
                mfc = 'r';
            else
                mk  = 'o';
                mfc = colors(i, :);
            end
            scatter(f1, f2, ptArea, 'Marker', mk, 'MarkerFaceColor', mfc, ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', strrep(expName, '_', '\_'));
            plottedNames{end+1} = expName; %#ok<AGROW>
        catch me
            warning('处理实验 %s 时出错: %s', expName, me.message);
        end
    end
    hold off;
    xlabel('f1', 'FontSize', 11, 'FontName', fontName);
    ylabel('f2', 'FontSize', 11, 'FontName', fontName);
    if ~isempty(plottedNames)
        leg = legend(plottedNames, 'Interpreter', 'none', 'Location', 'best');
        set(leg, 'FontName', fontName, 'FontSize', 10);
    end
    grid off;
    set(gca, 'FontSize', 10, 'FontName', fontName, 'Box', 'on', 'LineWidth', axisLineWidth);

    if ~isempty(plottedNames)
        figName = fullfile(ablation_dir, ['all_HVrank' num2str(hv_rank) '_obj_scatter']);
        saveas(gcf, [figName '.png']);
        savefig(gcf, [figName '.fig']);
        fprintf('已保存: %s.png / .fig（HV 排名第 %d，共 %d 个实验）\n', figName, hv_rank, length(plottedNames));
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

function runIdx = get_run_index_by_HV_rank_xlsx(summaryPath, hv_rank)
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

function Obj = read_obj_data_xlsx(runPath)
    if exist('readmatrix', 'file')
        Obj = readmatrix(runPath, 'Sheet', 'obj_data');
    else
        Obj = xlsread(runPath, 'obj_data');
    end
    if isempty(Obj) || size(Obj, 2) < 2
        Obj = [];
    end
end

function [runIdx, Obj] = get_obj_by_HV_rank_mat(matPath, hv_rank)
    S = load(matPath, 'runs');
    runs = S.runs;
    HV_end = arrayfun(@(s)s.hv_end, runs(:));
    [~, ord] = sort(HV_end, 'descend');
    idx = min(hv_rank, numel(ord));
    runIdx = runs(ord(idx)).Run;
    Obj = runs(ord(idx)).Obj_end;
    if isempty(Obj) || size(Obj,2) < 2
        Obj = [];
    end
end
