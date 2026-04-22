% =====================================================================
% 消融实验（并行版）：多次独立运行，保存为 .mat（不再写 Excel）
% 支持多算法依次运行，每个算法输出 1 份 mat：
%   Order_Data/ablation/<expName>_<算法名>/<expName>_<算法名>_runs.mat
%
% mat 文件结构：
%   runs(1..numRuns) : struct，每次运行一条，字段尽量对齐原 Excel 的 sheet：
%     - Obj_end      (N×M)
%     - Dec_end      (N×D)
%     - hv_curve     (nGen×1)
%     - igd_curve    (nGen×1)
%     - hv_end       (scalar)
%     - igd_end      (scalar)
%     - clustTag_end (N×1)  % 若算法未写 add_clustTag，则默认 1..N
%     - Details      (N×(2+2*m+1))  % pro.CalDetails(Dec_end) 输出
%     - run_time_sec (scalar) % 该次独立运行的墙钟时间（秒），parfor 下为各 worker 上单次耗时
%   summary : struct
%     - hv_end_all, igd_end_all, hv_curves, igd_curves, hv_mean, hv_std, obj_all_runs
%     - run_time_sec_all, run_time_mean, run_time_std, run_time_min, run_time_max
% 距离计算模式（在下方设置）：
%   global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = true;   % 含地图惩罚（默认）
%   global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = false;  % 纯欧氏距离
%
% NSGAIII_G 的 OperatorAMO 参数 pc / er（与 OperatorAMO.m 默认 0.5 / 0.25 一致）：
%   不必改 OperatorAMO.m；NSGAIII_G 内 [pc,er]=Algorithm.ParameterSet(0.5,0.25)，
%   本脚本用 platemo('algorithm', {@NSGAIII_G, pc, er}) 覆盖。
%   开启 doPcErSweep 时对 (pcList×erList) 网格各存一份 mat。
% =====================================================================

% --- 消融实验参数（修改此处）---
numRuns       = 31;    % 独立运行次数（如 20 或 30）
expName       = 'ablation0415';   % 实验基础名称，输出文件名会变为 <expName>_<算法名>
problem       = @planner_order_v3;
% 算法：本脚本仅做 NSGAIII_G 的 pc/er 参数实验（与 NSGAIII_G.m 中 ParameterSet 一致）
algorithm     = @NSGAIII_G;
N             = 300;   % 种群规模（与 platemo 一致）
maxFE         = 180000;
saveInterval  = 30;    % platemo 保存间隔
M             = 2;     % 目标数
plotLastRun   = false; % 是否在最后一轮画聚类散点图
% 输出根目录（可改为任意你想保存的路径；留空则用默认）
% 例：outRootDir = 'D:\AMOS_OUT\ablation';
outRootDir    = '/home/haichao/Documents/why/why_fromlev/PlatEMO/forOrderNew26/Order_Data/newdata2604/ablation0415';    % 默认：<本脚本目录>/../Order_Data/ablation
% --- 可选：切换距离模式 ---
% global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = true;

% --- OperatorAMO（pc / er）参数实验：网格 (pcList × erList) ---
doPcErSweep   = true;           % true：对每个 (pc,er) 组合各存一份 *_runs.mat
pcList        = [0.25, 0.5, 0.75];
erList        = [0.25, 0.5, 0.75];
defaultPc     = 0.5;             % doPcErSweep=false 时注入 NSGAIII_G
defaultEr     = 0.25;

% 统一为 cell，支持单算法或多算法依次运行
if ~iscell(algorithm), algorithm = {algorithm}; end

if doPcErSweep
    [PP, EE] = ndgrid(pcList(:), erList(:));
    pcErGrid = [PP(:), EE(:)];
else
    pcErGrid = [defaultPc, defaultEr];
end

for a = 1:length(algorithm)
    alg       = algorithm{a};
    algTag    = func2str(alg);   % 算法标识，用于输出文件名
    canInj    = algorithm_supports_operamo_pc_er(alg);
    if doPcErSweep && ~canInj
        nComb = 1;
    else
        nComb = size(pcErGrid, 1);
    end
    for g = 1:nComb
        if canInj && (doPcErSweep || nComb == 1)
            pcUse = pcErGrid(g, 1);
            erUse = pcErGrid(g, 2);
            algPlatemo = build_algorithm_platemo_spec(alg, pcUse, erUse);
            if doPcErSweep
                sufPcEr = sprintf('_pc%g_er%g', pcUse, erUse);
            else
                sufPcEr = '';
            end
        else
            pcUse = NaN;
            erUse = NaN;
            algPlatemo = alg;
            sufPcEr = '';
        end
        expNameCur = [expName '_' algTag sufPcEr];
        fprintf('\n>> 正在运行算法: %s%s（输出前缀: %s）\n', algTag, sufPcEr, expNameCur);
        [summaryTable, hv_end_all, igd_end_all] = ...
            run_ablation_exp_parallel_mat(problem, algPlatemo, numRuns, N, maxFE, saveInterval, M, expNameCur, outRootDir, pcUse, erUse);
        fprintf('\n===== %s 汇总（%d 次运行）=====\n', algTag, numRuns);
        disp(summaryTable);
        fprintf('HV_end:  %.6f ± %.6f\n', mean(hv_end_all), std(hv_end_all));
        fprintf('IGD_end: %.6f ± %.6f\n', mean(igd_end_all), std(igd_end_all));
        if ismember('Run_time_sec', summaryTable.Properties.VariableNames)
            rt = summaryTable.Run_time_sec;
            fprintf('单次运行耗时 Run_time_sec: 均值=%.2f s, 标准差=%.2f s, 最小=%.2f s, 最大=%.2f s\n', ...
                mean(rt), std(rt), min(rt), max(rt));
        end
        fprintf('====================================\n\n');
    end
end


function tf = algorithm_supports_operamo_pc_er(alg)
% 仅本脚本固定使用 NSGAIII_G（[pc,er]=ParameterSet(0.5,0.25)）。
    tf = strcmp(func2str(alg), 'NSGAIII_G');
end


function spec = build_algorithm_platemo_spec(alg, pc, er)
% NSGAIII_G：platemo('algorithm', {@NSGAIII_G, pc, er})
    spec = {alg, pc, er};
end


function [summaryTable, hv_end_all, igd_end_all] = ...
    run_ablation_exp_parallel_mat(problem, algorithm, numRuns, N, maxFE, saveInterval, M, expName, outRootDir, operamo_pc, operamo_er)
% 并行多次运行：输出 1 个 runs.mat（包含 runs(1..numRuns) 与 summary）。
% operamo_pc / operamo_er：若传入且非 NaN，写入 summary（用于参数实验追溯）。

    script_dir = fileparts(mfilename('fullpath'));
    if nargin < 9 || isempty(outRootDir)
        outRootDir = fullfile(script_dir, '..', 'Order_Data', 'ablation');
    end
    if nargin < 11 || isempty(operamo_pc), operamo_pc = NaN; end
    if nargin < 12 || isempty(operamo_er), operamo_er = NaN; end
    output_dir = fullfile(outRootDir, expName);
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    runs(numRuns,1) = struct('Run',[],'Obj_end',[],'Dec_end',[],'hv_curve',[],'igd_curve',[], ...
                             'hv_end',[],'igd_end',[],'clustTag_end',[],'Details',[],'run_time_sec',[]);

    % 为可复现：每次运行固定一个 seed（与并行调度无关）
    seedBase = 1000;

    fprintf('[消融-并行] 启动 %d 次运行...\n', numRuns);
    parfor r = 1:numRuns
        rng(seedBase + r, 'twister');
        tRun = tic;
        [Dec_end, Obj_end, igd, hv, igd_end, hv_end, clustTag_end, Details] = ...
            runplatemo_all1_single_min(problem, algorithm, N, maxFE, saveInterval, M);
        runs(r).Run        = r;
        runs(r).Obj_end    = Obj_end;
        runs(r).Dec_end    = Dec_end;
        runs(r).hv_curve   = hv(:);
        runs(r).igd_curve  = igd(:);
        runs(r).hv_end     = hv_end;
        runs(r).igd_end    = igd_end;
        runs(r).clustTag_end = clustTag_end;
        runs(r).Details    = Details;
        runs(r).run_time_sec = toc(tRun);
    end

    hv_end_all  = arrayfun(@(s)s.hv_end, runs(:));
    igd_end_all = arrayfun(@(s)s.igd_end, runs(:));
    run_time_sec_all = arrayfun(@(s)s.run_time_sec, runs(:));

    % 对齐曲线长度（列=run，行=代）
    maxGen = max(arrayfun(@(s)numel(s.hv_curve), runs(:)));
    hv_curves  = NaN(maxGen, numRuns);
    igd_curves = NaN(maxGen, numRuns);
    for r = 1:numRuns
        nGen = numel(runs(r).hv_curve);
        hv_curves(1:nGen, r)  = runs(r).hv_curve;
        igd_curves(1:nGen, r) = runs(r).igd_curve;
    end

    hv_mean = nanmean(hv_curves, 2);
    hv_std  = nanstd(hv_curves, 0, 2);

    obj_all_runs = [];
    for r = 1:numRuns
        Obj_end = runs(r).Obj_end;
        obj_all_runs = [obj_all_runs; [repmat(r, size(Obj_end,1), 1), Obj_end]]; %#ok<AGROW>
    end

    summary = struct();
    summary.expName      = expName;
    summary.operamo_pc   = operamo_pc;
    summary.operamo_er   = operamo_er;
    summary.outDir       = output_dir;
    summary.numRuns      = numRuns;
    summary.hv_end_all   = hv_end_all(:);
    summary.igd_end_all  = igd_end_all(:);
    summary.hv_curves    = hv_curves;
    summary.igd_curves   = igd_curves;
    summary.hv_mean      = hv_mean;
    summary.hv_std       = hv_std;
    summary.obj_all_runs = obj_all_runs;
    summary.run_time_sec_all = run_time_sec_all(:);
    summary.run_time_mean    = mean(run_time_sec_all);
    summary.run_time_std     = std(run_time_sec_all);
    summary.run_time_min     = min(run_time_sec_all);
    summary.run_time_max     = max(run_time_sec_all);

    summaryTable = table((1:numRuns)', hv_end_all(:), igd_end_all(:), run_time_sec_all(:), ...
        'VariableNames', {'Run', 'HV_end', 'IGD_end', 'Run_time_sec'});

    matPath = fullfile(output_dir, [expName '_runs.mat']);
    save(matPath, 'runs', 'summary', '-v7.3');
    fprintf('[消融-并行] 已保存: %s\n', matPath);
    fprintf('[消融-并行] 单次运行耗时: 均值=%.2f s, 标准差=%.2f s, 范围 [%.2f, %.2f] s\n', ...
        summary.run_time_mean, summary.run_time_std, summary.run_time_min, summary.run_time_max);
end


function [Dec_end, Obj_end, igd, hv, igd_end, hv_end, clustTag_end, Details] = ...
    runplatemo_all1_single_min(problem, algorithm, N, maxFE, saveInterval, M)
% 单次运行：调用 platemo，提取终点解与每代 HV/IGD（不写 Excel）。
% algorithm：函数句柄，或 cell（如 {@NSGAIII_G, 0.5, 0.25}），与 platemo 约定一致。

    platemo('problem', problem, 'algorithm', algorithm, 'N', N, 'M', M, 'maxFE', maxFE, 'save', saveInterval);
    result = evalin('base', 'result');
    pro = problem('M', M);

    lastSolution = result{end, end};
    num_obj = numel(lastSolution(1).obj);
    num_dec = numel(lastSolution(1).dec);
    Obj_end = zeros(numel(lastSolution), num_obj);
    Dec_end = zeros(numel(lastSolution), num_dec);
    clustTag_end = zeros(numel(lastSolution), 1);
    for idx = 1:numel(lastSolution)
        Obj_end(idx,:) = lastSolution(idx).obj;
        Dec_end(idx,:) = lastSolution(idx).dec;
        try
            ct = lastSolution(idx).add_clustTag;
            if ~isempty(ct), clustTag_end(idx) = ct(1); else, clustTag_end(idx) = idx; end
        catch
            clustTag_end(idx) = idx;
        end
    end

    igd_end = pro.CalMetric('IGD', result{end});
    hv_end  = pro.CalMetric('HV', result{end});

    nGen = size(result, 1);
    igd  = zeros(nGen, 1);
    hv   = zeros(nGen, 1);
    for i = 1:nGen
        igd(i) = pro.CalMetric('IGD', result{i, end});
        hv(i)  = pro.CalMetric('HV', result{i, end});
    end

    % Details：与原 Excel details_data 对齐
    Details = pro.CalDetails(Dec_end);
end


function [Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(problem,algorithm,N,maxFE,s,M)
platemo('problem',problem,'algorithm',algorithm,'N',N,'M',M,'maxFE',maxFE,'save',s);
result = evalin('base', 'result');

lastSolution = result{end, end};% 获取最后一行最后一列的 SOLUTION 结构
a=numel(lastSolution);
num_obj_columns = numel(lastSolution(1).obj);% 获取每个solution的obj值的列数
num_dec_columns = numel(lastSolution(1).dec);
% num_clustName_columns = numel(lastSolution(1).add_clustName);
% num_clustTag_columns = numel(lastSolution(1).add_clustTag);
% num_bj_columns = numel(lastSolution(1).bj);

Obj_end = zeros(numel(lastSolution),num_obj_columns); % 初始化存储obj值的数组
Dec_end = zeros(numel(lastSolution),num_dec_columns);
clustTag_end = zeros(numel(lastSolution), 1);  % 每个解的聚类编号
for idx = 1:numel(lastSolution)
    solution = lastSolution(idx);
    Obj = solution.obj;
    Obj_end(idx,:) = Obj;
    dec = solution.dec;
    Dec_end(idx,:) = dec;
    % 提取聚类标签（add_clustTag 字段，由 A_amos.m 写入）
    try
        ct = solution.add_clustTag;
        if ~isempty(ct)
            clustTag_end(idx) = ct(1);
        else
            clustTag_end(idx) = idx;
        end
    catch
        clustTag_end(idx) = idx;
    end
end
% num_pop_columns = numel(lastSolution(1).dec);
% Pop = zeros(numel(lastSolution),num_pop_columns); % 初始化存储obj值的数组。
% for idx = 1:numel(lastSolution)
%     solution = lastSolution(idx);
%     pop = solution.dec; % 提取每个1*1的solution的obj值
%     Pop(idx,:) = pop;
% end

% pro = problem('M',M,'dots',dots);
pro = problem('M',M);
igd_end=pro.CalMetric('IGD',result{end});
hv_end=pro.CalMetric('HV',result{end});
% 提取每代的IGD值
igd = zeros(size(result, 1), 1);
hv = zeros(size(result, 1), 1);
for i = 1:size(result, 1)
    igd(i) = pro.CalMetric('IGD', result{i, end});
    hv(i) = pro.CalMetric('HV', result{i, end});
end

objgen = zeros(numel(lastSolution), num_obj_columns);
Obj = cell(size(result, 1), 1);
for j = 1:size(result, 1)
for i = 1:numel(result{j, end})
    objgen = zeros(numel(result{j, end}), num_obj_columns);
    objgen(i,:) =result{j,2}(1,i).obj;    
end
    
    Obj{j,1}={objgen};
end



% popgen = zeros(numel(lastSolution), num_dec_columns);
Dec = cell(size(result, 1), 1);
for j = 1:size(result, 1)
for i = 1:numel(result{j, end})
    popgen = zeros(numel(result{j, end}), num_dec_columns);
    popgen(i,:) =result{j,2}(1,i).dec;    
end
    
    Dec{j,1}={popgen};
end

% clusttaggen = zeros(numel(lastSolution), num_clustTag_columns);
% clustTag = cell(size(result, 1), 1);
% for j = 2:size(result, 1)
% for i = 1:numel(result{j, end})
%     clusttaggen = zeros(numel(result{j, end}), num_clustTag_columns);
%     clusttaggen(i,:) =result{j,2}(1,i).add_clustTag;    
% end
% 
%     clustTag{j,1}={clusttaggen};
% end
% 
% clustTag = cell(size(result, 1), 1);
% for j = 1:size(result, 1)
% clusttaggen = zeros(numel(result{j, end}), num_clustTag_columns);
% for i = 1:numel(result{j, end})
%     clusttaggen(i,:) =result{j,2}(1,i).add_clustTag;    
% end   
%     clustTag{j,1}={clusttaggen};
% end
% 
% %clustnamegen = zeros(numel(lastSolution),num_clustName_columns);
% clustName = cell(size(result, 1), 1);
% for j = 1:size(result, 1)
%    clustnamegen = zeros(numel(result{j, end}),num_clustName_columns);
% for i = 1:numel(result{j, end})
%     clustnamegen(i,:) =result{j,2}(1,i).add_clustName;    
% end
%     clustName{j,1}={clustnamegen};
% end

% validtrait = cell(size(result, 1), 1);
% for j = 1:size(result, 1)
%     validtraitgen = zeros(numel(result{j, end}), 1);
% for i = 1:numel(result{j, end})
%     validtraitgen(i,:) =result{j,2}(1,i).validtrait;    
% end
% 
%     validtrait{j,1}={validtraitgen};
% end




% % 指定要保存的Excel文件名
% filename = 'maxhjddot2_disselect_initialselect_objduiying01.xlsx';
% 将 dot1 和 dot2 转换为字符串，用于文件名
% dot1_str = num2str(dot1);
% dot2_str = num2str(dot2);

% 指定保存目录为 Order_Data（forOrderNew26 下的 Order_Data 文件夹）
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', 'Order_Data','testData');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

% 创建唯一的文件名（保存到 Order_Data 中）
filename = fullfile(output_dir, '2026022224_09dma_order.xlsx');

% 将 obj 矩阵保存到 'Sheet1'
writematrix(Obj_end, filename, 'Sheet', 'obj_data');

% 将 dec 矩阵保存到 'Sheet2'
writematrix(Dec_end, filename, 'Sheet', 'dec_data');

% writematrix(validtrait_end, filename, 'Sheet', 'validtrait_data');
% 
% writematrix(bj_end, filename, 'Sheet', 'bj_data');

writematrix(hv, filename, 'Sheet', 'hv_data');

% 保存聚类标签（每行：f1, f2, clustTag）
clustData = [Obj_end, clustTag_end];
clustHeader = {'f1_energy_J', 'f2_time_s', 'clustTag'};
writecell(clustHeader, filename, 'Sheet', 'clust_data', 'Range', 'A1');
writematrix(clustData, filename, 'Sheet', 'clust_data', 'Range', 'A2');

% ---- 诊断：调用 CalDetails 获取每架无人机能耗/用时/惩罚项 ----
m = pro.m;  % 无人机数量
Details = pro.CalDetails(Dec_end);
% Details 列含义（m=3时）：
%   col 1       : Penalty1（载重超限，J当量）
%   col 2       : Penalty2（电量超限，J）
%   col 3       : E_drone1 (J)
%   col 4       : E_drone2 (J)
%   col 5       : E_drone3 (J)
%   col 6       : T_drone1 (s)
%   col 7       : T_drone2 (s)
%   col 8       : T_drone3 (s)
%   col 9       : 是否可行（1=可行，0=违约）

% 构建带表头的 cell 写入 Excel
colNames = ['Penalty1_J', 'Penalty2_J', ...
    arrayfun(@(k) sprintf('E_drone%d_J',k), 1:m, 'UniformOutput',false), ...
    arrayfun(@(k) sprintf('T_drone%d_s',k), 1:m, 'UniformOutput',false), ...
    'Feasible'];
header = [{'f1_energy_J','f2_time_s'}, colNames];
dataWithObj = [num2cell(Obj_end), num2cell(Details)];
detailSheet = [header; dataWithObj];
writecell(detailSheet, filename, 'Sheet', 'details_data');

% 控制台打印摘要
nFeasible = sum(Details(:, end));
fprintf('\n===== Pareto 前沿解诊断（共 %d 个解）=====\n', size(Dec_end,1));
fprintf('  可行解数量：%d / %d（%.1f%%）\n', nFeasible, size(Dec_end,1), 100*nFeasible/size(Dec_end,1));
fprintf('  载重违约（Penalty1>0）：%d 个解\n', sum(Details(:,1)>0));
fprintf('  电量违约（Penalty2>0）：%d 个解\n', sum(Details(:,2)>0));
fprintf('\n  可行解能耗范围：%.1f ~ %.1f kJ\n', ...
    min(sum(Details(Details(:,end)==1, 3:2+m), 2))/1000, ...
    max(sum(Details(Details(:,end)==1, 3:2+m), 2))/1000);
fprintf('  各无人机最大单机能耗（全部解）：');
for k = 1:m
    fprintf('Drone%d: %.1f kJ  ', k, max(Details(:,2+k))/1000);
end
fprintf('\n');
fprintf('  maxEC=800 kJ，超出该值即违约。\n');
fprintf('==========================================\n\n');

% 绘制聚类散点图
plot_clust_scatter(Obj_end, clustTag_end, get_filename_only_local(filename));
end

function name = get_filename_only_local(path)
    [~, name, ext] = fileparts(path);
    name = [name ext];
end


