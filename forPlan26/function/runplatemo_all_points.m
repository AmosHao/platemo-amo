%platemo('problem',@IMMOEA,'algorithm',@GLT1,'D',3,'N',100,'maxFE',10000,'save',20);
% clear
% [validtrait_end,validtrait,clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(@planner_points,@GOCEAfortravel_singlepath_3,100,20000,20,3,[0.2,0.2,0.8,3.5,3,0.8]);
% function [validtrait_end,validtrait,clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(problem,algorithm,N,maxFE,s,M,parameterdot1_x,dot1_y,dot1_z,dot2_x,dot2_y,dot2_z)
% platemo('problem',problem,'algorithm',algorithm,'N',N,'M',M,'maxFE',maxFE,'save',s,'dot1_x',dot1_x,'dot1_y',dot1_y,'dot1_z',dot1_z,'dot2_x',dot2_x,'dot2_y',dot2_y,'dot2_z',dot2_z);
function [validtrait_end,validtrait,clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all_points(problem,algorithm,N,maxFE,s,M,dots)
platemo('problem',problem,'algorithm',algorithm,'N',N,'M',M,'maxFE',maxFE,'save',s,'dots',dots);
result = evalin('base', 'result');

lastSolution = result{end, end};% 获取最后一行最后一列的 SOLUTION 结构
a=numel(lastSolution);
num_obj_columns = numel(lastSolution(1).obj);% 获取每个solution的obj值的列数
num_dec_columns = numel(lastSolution(1).dec);
num_clustName_columns = numel(lastSolution(1).add_clustName);
num_clustTag_columns = numel(lastSolution(1).add_clustTag);

Obj_end = zeros(numel(lastSolution),num_obj_columns); % 初始化存储obj值的数组
Dec_end = zeros(numel(lastSolution),num_dec_columns);
clustTag_end = zeros(numel(lastSolution),num_clustTag_columns);
clustName_end = zeros(numel(lastSolution),num_clustName_columns);
validtrait_end = zeros(numel(lastSolution),1);
for idx = 1:numel(lastSolution)
    solution = lastSolution(idx);
    Obj = solution.obj; % 提取每个1*1的solution的obj值
    Obj_end(idx,:) = Obj;
    dec = solution.dec;
    Dec_end(idx,:) = dec;
    add_clustTag = solution.add_clustTag;
    clustTag_end(idx,:) = add_clustTag;
    add_clustName = solution.add_clustName;
    clustName_end(idx,:) = add_clustName;
    validtrait = solution.validtrait;
    validtrait_end(idx,:) = validtrait;
end
% num_pop_columns = numel(lastSolution(1).dec);
% Pop = zeros(numel(lastSolution),num_pop_columns); % 初始化存储obj值的数组。
% for idx = 1:numel(lastSolution)
%     solution = lastSolution(idx);
%     pop = solution.dec; % 提取每个1*1的solution的obj值
%     Pop(idx,:) = pop;
% end

pro = problem('M',M);
% pro = problem{1,1};
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

clustTag = cell(size(result, 1), 1);
for j = 1:size(result, 1)
clusttaggen = zeros(numel(result{j, end}), num_clustTag_columns);
for i = 1:numel(result{j, end})
    clusttaggen(i,:) =result{j,2}(1,i).add_clustTag;    
end   
    clustTag{j,1}={clusttaggen};
end

%clustnamegen = zeros(numel(lastSolution),num_clustName_columns);
clustName = cell(size(result, 1), 1);
for j = 1:size(result, 1)
   clustnamegen = zeros(numel(result{j, end}),num_clustName_columns);
for i = 1:numel(result{j, end})
    clustnamegen(i,:) =result{j,2}(1,i).add_clustName;    
end
    clustName{j,1}={clustnamegen};
end

validtrait = cell(size(result, 1), 1);
for j = 1:size(result, 1)
    validtraitgen = zeros(numel(result{j, end}), 1);
for i = 1:numel(result{j, end})
    validtraitgen(i,:) =result{j,2}(1,i).validtrait;    
end
    
    validtrait{j,1}={validtraitgen};
end
% 
% % 指定要保存的Excel文件名
% filename = 'point2_penal2_dmin0.3_travel_singlepath1015_F3_penaly0_guideclose_maxfe20000_initialy_rand0.3.xlsx';
% 
% % 将 obj 矩阵保存到 'Sheet1'
% writematrix(Obj_end, filename, 'Sheet', 'obj_data');
% 
% % 将 dec 矩阵保存到 'Sheet2'
% writematrix(Dec_end, filename, 'Sheet', 'dec_data');
% 
% writematrix(validtrait_end, filename, 'Sheet', 'validtrait_data');
end


