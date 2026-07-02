% =====================================================================
% 消融实验：多次独立运行，每轮保存一份完整 Excel，并生成汇总文件便于画图
% 支持多算法依次运行，输出文件名带算法标识：<expName>_<算法名>_run01.xlsx 等
% =====================================================================
% 输出结构（目录 Order_Data/ablation/<expName>_<算法名>/）：
%   1) 每轮一份完整 Excel：<expName>_<算法名>_run01.xlsx ~ run30.xlsx
%      与单次运行结构相同：obj_data, dec_data, hv_data, igd_data, clust_data, details_data
%   2) 汇总 <expName>_<算法名>_summary.xlsx：
%      - runs：每行一次运行，列 Run, HV_end, IGD_end
%      - stats：HV_end/IGD_end 的 Mean, Std
%      - hv_curves：行=代数、列=第几次运行（读作 30 条 HV 曲线）
%      - igd_curves：同上
%      - hv_mean_curve：表头 Gen, HV_mean, HV_std → 直接画 30 次平均 HV 曲线（可加 ±std 阴影）
%      - obj_all_runs：表头 Run, f1, f2，行为 30*N 条终点解 → 画总前沿面（可先取非支配再画）
% 距离计算模式（在下方设置）：
%   global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = true;   % 含地图惩罚（默认）
%   global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = false;  % 纯欧氏距离
% =====================================================================

% --- 消融实验参数（修改此处）---
numRuns       = 1;    % 独立运行次数（如 20 或 30）
expName       = 'fangzhen7';   % 实验基础名称，输出文件名会变为 <expName>_<算法名>
problem       = @planner_order_v3_fangzhen;
% 算法：可写单个句柄，或 cell 数组依次运行多个（输出带算法标识）
% algorithm     = @A_amos;
% algorithm     = {@A_amos, @A_amos_noCluster, @A_amos_allBYJC};  % 多算法示例
algorithm     = {@NSGAIII_G};
N             = 300;   % 种群规模（与 platemo 一致）
maxFE         = 180000;
saveInterval  = 30;    % platemo 保存间隔
M             = 2;     % 目标数
plotLastRun   = false; % 是否在最后一轮画聚类散点图
% --- 可选：切换距离模式 ---
% global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = false;

% 统一为 cell，支持单算法或多算法依次运行
if ~iscell(algorithm), algorithm = {algorithm}; end

for a = 1:length(algorithm)
    alg       = algorithm{a};
    algTag    = func2str(alg);   % 算法标识，用于输出文件名
    expNameCur = [expName '_' algTag];  % 带算法标识的实验名，如 ablation_A_amos
    fprintf('\n>> 正在运行算法: %s（输出前缀: %s）\n', algTag, expNameCur);
    [summaryTable, hv_end_all, igd_end_all, hv_curves, igd_curves] = ...
        run_ablation_exp(problem, alg, numRuns, N, maxFE, saveInterval, M, expNameCur, plotLastRun);
    fprintf('\n===== %s 汇总（%d 次运行）=====\n', algTag, numRuns);
    disp(summaryTable);
    fprintf('HV_end:  %.6f ± %.6f\n', mean(hv_end_all), std(hv_end_all));
    fprintf('IGD_end: %.6f ± %.6f\n', mean(igd_end_all), std(igd_end_all));
    fprintf('====================================\n\n');
end


function [summaryTable, hv_end_all, igd_end_all, hv_curves, igd_curves] = ...
    run_ablation_exp(problem, algorithm, numRuns, N, maxFE, saveInterval, M, expName, plotLastRun)
% 多次独立运行：每轮保存一份完整 Excel（run01~run30，与单次运行结构一致），并生成汇总文件。
% 汇总文件含：runs/stats、hv_curves/igd_curves（每列一次运行）、hv_mean_curve（画平均HV曲线）、
% obj_all_runs（30次全部解，用于画总前沿面）。
    hv_end_all   = zeros(numRuns, 1);
    igd_end_all  = zeros(numRuns, 1);
    hv_curves    = [];
    igd_curves   = [];
    obj_all_runs = [];  % [Run, f1, f2, ...] 所有运行的终点目标合并，用于总前沿面

    script_dir   = fileparts(mfilename('fullpath'));
    output_dir   = fullfile(script_dir, '..', 'Order_Data', 'ablation', expName);
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    for r = 1:numRuns
        fprintf('[消融] 第 %d/%d 次运行...\n', r, numRuns);
        [Dec, Obj, Dec_end, Obj_end, igd, hv, igd_end, hv_end, result, pro] = ...
            runplatemo_all1_single(problem, algorithm, N, maxFE, saveInterval, M);

        hv_end_all(r)  = hv_end;
        igd_end_all(r) = igd_end;
        % 合并终点目标：每行 [Run, f1, f2]
        obj_all_runs = [obj_all_runs; [repmat(r, size(Obj_end,1), 1), Obj_end]];

        nGen = length(hv);
        if isempty(hv_curves)
            hv_curves  = hv(:);
            igd_curves = igd(:);
        else
            nOld = size(hv_curves, 1);
            if nGen > nOld
                hv_curves  = [hv_curves; NaN(nGen - nOld, size(hv_curves, 2))];
                igd_curves = [igd_curves; NaN(nGen - nOld, size(igd_curves, 2))];
            end
            hv_curves(1:nGen, r)  = hv(1:nGen);
            igd_curves(1:nGen, r) = igd(1:nGen);
            if nGen < size(hv_curves, 1)
                hv_curves(nGen+1:end, r)  = NaN;
                igd_curves(nGen+1:end, r) = NaN;
            end
        end

        % 每轮都保存一份完整 Excel（与 runplatemo_amo 单次运行结构一致）
        filename = fullfile(output_dir, [expName '_run' sprintf('%02d', r) '.xlsx']);
        save_single_run_results(Dec_end, Obj_end, hv, igd, result, pro, filename, ...
            [expName '_run' sprintf('%02d', r)], r == numRuns && plotLastRun);
    end

    if size(hv_curves, 2) < numRuns
        hv_curves(:, end+1:numRuns)  = NaN;
        igd_curves(:, end+1:numRuns) = NaN;
    end

    % 汇总表
    summaryTable = table((1:numRuns)', hv_end_all, igd_end_all, ...
        'VariableNames', {'Run', 'HV_end', 'IGD_end'});

    summaryPath = fullfile(output_dir, [expName '_summary.xlsx']);
    writetable(summaryTable, summaryPath, 'Sheet', 'runs');
    stats = table({'HV_end'; 'IGD_end'}, ...
        [mean(hv_end_all); mean(igd_end_all)], ...
        [std(hv_end_all); std(igd_end_all)], ...
        'VariableNames', {'Metric', 'Mean', 'Std'});
    writetable(stats, summaryPath, 'Sheet', 'stats');
    writematrix(hv_curves, summaryPath, 'Sheet', 'hv_curves');   % 列=run，行=代
    writematrix(igd_curves, summaryPath, 'Sheet', 'igd_curves');

    % 便于画“30次平均HV曲线”：Gen, HV_mean, HV_std
    nGen = size(hv_curves, 1);
    hv_mean = nanmean(hv_curves, 2);
    hv_std  = nanstd(hv_curves, 0, 2);
    hv_mean_curve = [(1:nGen)', hv_mean, hv_std];
    writecell({'Gen', 'HV_mean', 'HV_std'}, summaryPath, 'Sheet', 'hv_mean_curve', 'Range', 'A1');
    writematrix(hv_mean_curve, summaryPath, 'Sheet', 'hv_mean_curve', 'Range', 'A2');

    % 便于画“总前沿面”：所有运行的终点目标合并，列 Run, f1, f2
    writecell({'Run', 'f1', 'f2'}, summaryPath, 'Sheet', 'obj_all_runs', 'Range', 'A1');
    writematrix(obj_all_runs, summaryPath, 'Sheet', 'obj_all_runs', 'Range', 'A2');
end


function save_single_run_results(Dec_end, Obj_end, hv, igd, result, pro, filename, expName, doPlot)
    numSol = size(Obj_end, 1);
    clustTag_end = zeros(numSol, 1);
    lastSolution = result{end, end};
    for idx = 1:numel(lastSolution)
        try
            ct = lastSolution(idx).add_clustTag;
            if ~isempty(ct), clustTag_end(idx) = ct(1); else, clustTag_end(idx) = idx; end
        catch, clustTag_end(idx) = idx; end
    end
    writematrix(Obj_end, filename, 'Sheet', 'obj_data');
    writematrix(Dec_end, filename, 'Sheet', 'dec_data');
    writematrix(hv, filename, 'Sheet', 'hv_data');
    writematrix(igd, filename, 'Sheet', 'igd_data');
    clustData = [Obj_end, clustTag_end];
    writecell({'f1_energy_J','f2_time_s','clustTag'}, filename, 'Sheet', 'clust_data', 'Range', 'A1');
    writematrix(clustData, filename, 'Sheet', 'clust_data', 'Range', 'A2');
    m = pro.m;
    Details = pro.CalDetails(Dec_end);
    colNames = ['Penalty1_J','Penalty2_J', ...
        arrayfun(@(k) sprintf('E_drone%d_J',k), 1:m, 'UniformOutput',false), ...
        arrayfun(@(k) sprintf('T_drone%d_s',k), 1:m, 'UniformOutput',false), 'Feasible'];
    header = [{'f1_energy_J','f2_time_s'}, colNames];
    writecell([header; num2cell([Obj_end, Details])], filename, 'Sheet', 'details_data');
    if doPlot && exist('plot_clust_scatter','file')
        plot_clust_scatter(Obj_end, clustTag_end, [expName '_run_last']);
    end
end


function [Dec, Obj, Dec_end, Obj_end, igd, hv, igd_end, hv_end, result, pro] = ...
    runplatemo_all1_single(problem, algorithm, N, maxFE, saveInterval, M)
% 单次运行：调用 platemo，提取 Dec/Obj/每代及终点 IGD、HV。
    platemo('problem', problem, 'algorithm', algorithm, 'N', N, 'M', M, 'maxFE', maxFE, 'save', saveInterval);
    result = evalin('base', 'result');
    pro = problem('M', M);
    lastSolution = result{end, end};
    num_obj = numel(lastSolution(1).obj);
    num_dec = numel(lastSolution(1).dec);
    Obj_end = zeros(numel(lastSolution), num_obj);
    Dec_end = zeros(numel(lastSolution), num_dec);
    for idx = 1:numel(lastSolution)
        Obj_end(idx,:) = lastSolution(idx).obj;
        Dec_end(idx,:) = lastSolution(idx).dec;
    end
    igd_end = pro.CalMetric('IGD', result{end});
    hv_end  = pro.CalMetric('HV', result{end});
    nGen = size(result, 1);
    igd = zeros(nGen, 1);
    hv  = zeros(nGen, 1);
    for i = 1:nGen
        igd(i) = pro.CalMetric('IGD', result{i, end});
        hv(i)  = pro.CalMetric('HV', result{i, end});
    end
    Obj = cell(nGen, 1);
    Dec = cell(nGen, 1);
    for j = 1:nGen
        objgen = zeros(numel(result{j,end}), num_obj);
        popgen = zeros(numel(result{j,end}), num_dec);
        for i = 1:numel(result{j,end})
            objgen(i,:) = result{j,2}(1,i).obj;
            popgen(i,:) = result{j,2}(1,i).dec;
        end
        Obj{j} = {objgen};
        Dec{j} = {popgen};
    end
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


