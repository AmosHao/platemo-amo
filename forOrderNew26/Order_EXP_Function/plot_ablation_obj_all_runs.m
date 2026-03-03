% =====================================================================
% 绘制指定实验的 31 次运行全部目标值 (f1-f2) 散点图（不按 HV 排名筛选）
% 数据：Order_Data/ablation/<实验名>/ 下 <实验名>_run01.xlsx ~ run30/31 的 obj_data
% =====================================================================
% 用法：在编辑器中按 F5 运行；或命令窗口：
%   plot_ablation_obj_all_runs(expNames)
% 参数：
%   expNames - 实验名（与 ablation 下子文件夹一致），cell 或字符。留空则自动扫描所有实验
% 输出：Order_Data/ablation/<实验名>_all_runs_obj_scatter.png（及 .fig），该实验所有 run 在一张图
%       若指定多个实验，则输出 all_exps_all_runs_obj_scatter.png，每个实验一种颜色
% =====================================================================

% ---------- 直接运行本文件时执行（F5）：可修改下面一行 ----------
expNames = {'ablation06_A_amos', 'ablation05_A_amos_noCluster', 'ablation05_A_amos_allBYJC'};   % 留空=扫描所有实验
plot_ablation_obj_all_runs_inside(expNames);

function plot_ablation_obj_all_runs_inside(expNames)
    run_plot_ablation_obj_all_runs(expNames);
end

function run_plot_ablation_obj_all_runs(expNames)
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
    end

    % 每个实验一种颜色，一张图上画所有 run
    colors = lines(max(length(expNames), 1));

    figure('Name', '各实验全部 run 目标值', 'NumberTitle', 'off');
    hold on;
    plottedNames = {};
    for i = 1:length(expNames)
        expName = expNames{i};
        expPath = fullfile(ablation_dir, expName);
        if ~isfolder(expPath)
            warning('未找到实验目录: %s，跳过。', expPath);
            continue;
        end
        % 扫描 run01 ~ run31
        all_f1 = [];
        all_f2 = [];
        runCount = 0;
        for runIdx = 1:31
            runPath = fullfile(expPath, [expName '_run' sprintf('%02d', runIdx) '.xlsx']);
            if ~isfile(runPath)
                continue;
            end
            try
                Obj = read_obj_data(runPath);
                if isempty(Obj)
                    continue;
                end
                all_f1 = [all_f1; Obj(:, 1)]; %#ok<AGROW>
                all_f2 = [all_f2; Obj(:, 2)]; %#ok<AGROW>
                runCount = runCount + 1;
            catch
                % 单次 run 读取出错则跳过该 run
            end
        end
        if isempty(all_f1)
            warning('实验 %s 下未找到有效 run 的 obj_data，跳过。', expName);
            continue;
        end
        scatter(all_f1, all_f2, 18, colors(i,:), 'filled', 'DisplayName', sprintf('%s (%d runs)', strrep(expName, '_', '\_'), runCount));
        plottedNames{end+1} = expName; %#ok<AGROW>
    end
    hold off;
    xlabel('f1 (能耗)', 'FontSize', 11);
    ylabel('f2 (时间)', 'FontSize', 11);
    title('指定实验 — 全部 run 目标值（不按 HV 排名）');
    if ~isempty(plottedNames)
        legend(plottedNames, 'Interpreter', 'none', 'Location', 'best');
    end
    grid on;
    set(gca, 'FontSize', 10);

    if ~isempty(plottedNames)
        if length(plottedNames) == 1
            figName = fullfile(ablation_dir, plottedNames{1}, [plottedNames{1} '_all_runs_obj_scatter']);
        else
            figName = fullfile(ablation_dir, 'all_exps_all_runs_obj_scatter');
        end
        saveas(gcf, [figName '.png']);
        savefig(gcf, [figName '.fig']);
        fprintf('已保存: %s.png / .fig（共 %d 个实验，每实验所有 run）\n', figName, length(plottedNames));
    end
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
