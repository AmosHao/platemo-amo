% =====================================================================
% 分析多组消融实验的 summary 数据，对比 HV_end/IGD_end 与收敛曲线，辅助排查“完整版不如消融版”的原因
% 数据来源：Order_Data/ablation/<实验名>/<实验名>_summary.xlsx
% =====================================================================
% 用法：在编辑器中按 F5 运行；或命令窗口 analyze_ablation_summary(expNames)
% 参数 expNames：cell，如 {'ablation01_A_amos','ablation01_A_amos_noCluster','ablation01_A_amos_noAdapt'}
% =====================================================================

expNames = {'ablation06_A_amos', 'ablation06_A_amos_noCluster','ablation06_A_amos_allBYJC'};
analyze_ablation_summary_inside(expNames);

function analyze_ablation_summary_inside(expNames)
    if nargin < 1 || isempty(expNames)
        expNames = {'ablation01_A_amos', 'ablation01_A_amos_noCluster', 'ablation01_A_amos_noAdapt'};
    end
    if ischar(expNames), expNames = cellstr(expNames); end

    script_dir   = fileparts(mfilename('fullpath'));
    ablation_dir = fullfile(script_dir, '..', 'Order_Data', 'ablation');

    if ~isfolder(ablation_dir)
        fprintf('未找到目录: %s\n', ablation_dir);
        return;
    end

    nExp = length(expNames);
    HV_end_mean = zeros(nExp, 1);
    HV_end_std  = zeros(nExp, 1);
    HV_end_min  = zeros(nExp, 1);
    HV_end_max  = zeros(nExp, 1);
    IGD_end_mean = zeros(nExp, 1);
    IGD_end_std  = zeros(nExp, 1);
    lastGen_HV_mean = zeros(nExp, 1);
    lastGen_HV_std  = zeros(nExp, 1);
    hv_curves_all  = cell(nExp, 1);
    nGen_all       = zeros(nExp, 1);

    for i = 1:nExp
        expName = expNames{i};
        summaryPath = fullfile(ablation_dir, expName, [expName '_summary.xlsx']);
        if ~isfile(summaryPath)
            fprintf('未找到: %s\n', summaryPath);
            HV_end_mean(i) = NaN; HV_end_std(i) = NaN; IGD_end_mean(i) = NaN; IGD_end_std(i) = NaN;
            lastGen_HV_mean(i) = NaN; lastGen_HV_std(i) = NaN;
            continue;
        end

        % runs: Run, HV_end, IGD_end
        try
            T = readtable(summaryPath, 'Sheet', 'runs', 'VariableNamingRule', 'preserve');
        catch
            T = readtable(summaryPath, 'Sheet', 'runs');
        end
        vars = T.Properties.VariableNames;
        hvCol  = find(cellfun(@(x) ~isempty(regexpi(x, 'HV_end|HVend')), vars), 1);
        igdCol = find(cellfun(@(x) ~isempty(regexpi(x, 'IGD_end|IGDend')), vars), 1);
        if isempty(hvCol),  hvCol = 2; end
        if isempty(igdCol), igdCol = 3; end
        HV_end  = T.(vars{hvCol});
        IGD_end = T.(vars{igdCol});

        HV_end_mean(i) = mean(HV_end);
        HV_end_std(i)  = std(HV_end);
        HV_end_min(i)  = min(HV_end);
        HV_end_max(i)  = max(HV_end);
        IGD_end_mean(i) = mean(IGD_end);
        IGD_end_std(i)  = std(IGD_end);

        % hv_mean_curve: 最后一行的 HV_mean, HV_std
        try
            Tc = readtable(summaryPath, 'Sheet', 'hv_mean_curve', 'VariableNamingRule', 'preserve');
        catch
            Tc = readtable(summaryPath, 'Sheet', 'hv_mean_curve');
        end
        vc = Tc.Properties.VariableNames;
        mCol = find(cellfun(@(x) ~isempty(regexpi(x, 'HV_mean|HVmean')), vc), 1);
        sCol = find(cellfun(@(x) ~isempty(regexpi(x, 'HV_std|HVstd')), vc), 1);
        if isempty(mCol), mCol = 2; end
        if isempty(sCol), sCol = 3; end
        lastGen_HV_mean(i) = Tc.(vc{mCol})(end);
        lastGen_HV_std(i)  = Tc.(vc{sCol})(end);

        % hv_curves: 每列一次运行，行=代 → 用于看每代 30 次的方差
        try
            M = readmatrix(summaryPath, 'Sheet', 'hv_curves');
        catch
            M = xlsread(summaryPath, 'hv_curves');
        end
        hv_curves_all{i} = M;
        nGen_all(i) = size(M, 1);
    end

    %% 打印对比报告
    fprintf('\n========== 消融实验 summary 对比 ==========\n');
    fprintf('数据目录: %s\n\n', ablation_dir);

    fprintf('--- 1) 30 次运行终点 HV_end（越大越好）---\n');
    fprintf('%-35s  Mean      Std       Min       Max\n', '实验名');
    for i = 1:nExp
        fprintf('%-35s  %.6f  %.6f  %.6f  %.6f\n', ...
            expNames{i}, HV_end_mean(i), HV_end_std(i), HV_end_min(i), HV_end_max(i));
    end

    fprintf('\n--- 2) 30 次运行终点 IGD_end（越小越好）---\n');
    fprintf('%-35s  Mean      Std\n', '实验名');
    for i = 1:nExp
        fprintf('%-35s  %.6f  %.6f\n', expNames{i}, IGD_end_mean(i), IGD_end_std(i));
    end

    fprintf('\n--- 3) 最后一代 30 次平均 HV（与 hv_mean_curve 一致）---\n');
    fprintf('%-35s  HV_mean   HV_std(30次)\n', '实验名');
    for i = 1:nExp
        fprintf('%-35s  %.6f  %.6f\n', expNames{i}, lastGen_HV_mean(i), lastGen_HV_std(i));
    end

    % 每代 30 次 HV 的方差（标准差）— 看谁波动大
    fprintf('\n--- 4) 收敛过程稳定性：每代 30 次 HV 的标准差（代 1 / 中间 / 末代）---\n');
    fprintf('%-35s  Gen1_std  Mid_std   LastGen_std\n', '实验名');
    for i = 1:nExp
        M = hv_curves_all{i};
        if isempty(M) || size(M, 2) < 2
            fprintf('%-35s  (无数据)\n', expNames{i});
            continue;
        end
        perGenStd = nanstd(M, 0, 2);
        g1 = perGenStd(1);
        gMid = perGenStd(max(1, round(end/2)));
        gLast = perGenStd(end);
        fprintf('%-35s  %.6f  %.6f  %.6f\n', expNames{i}, g1, gMid, gLast);
    end

    % 简要诊断
    fprintf('\n--- 5) 简要诊断 ---\n');
    [~, bestIdx] = max(HV_end_mean);
    [~, worstIdx] = min(HV_end_mean);
    fprintf('平均 HV_end 最高: %s (%.6f)\n', expNames{bestIdx}, HV_end_mean(bestIdx));
    fprintf('平均 HV_end 最低: %s (%.6f)\n', expNames{worstIdx}, HV_end_mean(worstIdx));
    if HV_end_std(worstIdx) > HV_end_std(bestIdx) * 1.2
        fprintf('→ 完整版 30 次运行方差更大，表现不稳定，易受聚类/自适应影响。\n');
    end
    fprintf('============================================\n\n');
end
