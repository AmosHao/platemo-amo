%platemo('problem',@IMMOEA,'algorithm',@GLT1,'D',3,'N',100,'maxFE',10000,'save',20);
%clear

% =====================================================================
% 距离计算模式切换（在调用 runplatemo_all1 之前设置）：
%
%   模式A：含地图惩罚（禁飞区/人流密集区/噪音区），默认
%     global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = true;
%     （或不设置，类默认值即为 true）
%
%   模式B：纯欧氏距离（忽略所有障碍物惩罚，便于对比）
%     global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = false;
%
%   对比完毕后清除全局变量，恢复默认：
%     clear global AMOS_USE_MAP_PENALTY;
% =====================================================================

% --- 在此处取消注释以切换模式 ---
% global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = false;  % 纯欧氏距离
% global AMOS_USE_MAP_PENALTY; AMOS_USE_MAP_PENALTY = true;   % 含惩罚（默认）

[Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(@planner_order_v3,@A_amos,100,80000,20,2);
% [clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(@planner_travel_maxhjd02,@RMDE,100,10000,20,3,dots,dot1,dot2);
% end
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
% clustTag_end = zeros(numel(lastSolution),num_clustTag_columns);
% clustName_end = zeros(numel(lastSolution),num_clustName_columns);
% validtrait_end = zeros(numel(lastSolution),1);
% bj_end = zeros(numel(lastSolution),num_bj_columns);
for idx = 1:numel(lastSolution)
    solution = lastSolution(idx);
    Obj = solution.obj; % 提取每个1*1的solution的obj值
    Obj_end(idx,:) = Obj;
    dec = solution.dec;
    Dec_end(idx,:) = dec;
    % add_clustTag = solution.add_clustTag;
    % clustTag_end(idx,:) = add_clustTag;
    % add_clustName = solution.add_clustName;
    % clustName_end(idx,:) = add_clustName;
    % validtrait = solution.validtrait;
    % validtrait_end(idx,:) = validtrait;
    % bj = solution.bj;
    % bj_end(idx,:) = bj;
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
filename = fullfile(output_dir, '2026022221_dma_order.xlsx');

% 将 obj 矩阵保存到 'Sheet1'
writematrix(Obj_end, filename, 'Sheet', 'obj_data');

% 将 dec 矩阵保存到 'Sheet2'
writematrix(Dec_end, filename, 'Sheet', 'dec_data');

% writematrix(validtrait_end, filename, 'Sheet', 'validtrait_data');
% 
% writematrix(bj_end, filename, 'Sheet', 'bj_data');

writematrix(hv, filename, 'Sheet', 'hv_data');

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
end


