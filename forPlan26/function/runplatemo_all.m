%platemo('problem',@test1GOCEA_new,'algorithm',@GLT1,'D',3,'N',100,'maxFE',10000,'save',20);
%[clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(@planner_original,@PreDEMO,100,30000,20,30,3);
function [clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all(problem,algorithm,N,maxFE,s,vardim,M)
platemo('problem',problem,'algorithm',algorithm,'D',vardim,'N',N,'M',M,'maxFE',maxFE,'save',s);
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
end
% num_pop_columns = numel(lastSolution(1).dec);
% Pop = zeros(numel(lastSolution),num_pop_columns); % 初始化存储obj值的数组。
% for idx = 1:numel(lastSolution)
%     solution = lastSolution(idx);
%     pop = solution.dec; % 提取每个1*1的solution的obj值
%     Pop(idx,:) = pop;
% end

pro = problem('M',M,'D',vardim);
igd_end=pro.CalMetric('IGD',result{end});
hv_end=pro.CalMetric('HV',result{end});
% 提取每代的IGD值
igd = zeros(size(result, 1), 1);
hv = zeros(size(result, 1), 1);
for i = 1:size(result, 1)
    igd(i) = pro.CalMetric('IGD', result{i, end});
    hv(i) = pro.CalMetric('HV', result{i, end});
end



Obj = cell(size(result, 1), 1);
for j = 1:size(result, 1)
    genSolution = result{j, 2};
    objgen = zeros(numel(genSolution), num_obj_columns);
for i = 1:numel(genSolution)
    objgen(i,:) =result{j,2}(1,i).obj;    
end
    
    Obj{j,1}={objgen};
end




Dec = cell(size(result, 1), 1);
for j = 1:size(result, 1)
    genSolution = result{j, 2};
    popgen = zeros(numel(genSolution), num_dec_columns);

for i = 1:numel(genSolution)
    popgen(i,:) =result{j,2}(1,i).dec;    
end
    
    Dec{j,1}={popgen};
end


clustTag = cell(size(result, 1), 1);
for j = 1:size(result, 1)
    genSolution = result{j, 2};
    clusttaggen = zeros(numel(genSolution), num_clustTag_columns);
for i = 1:numel(genSolution)
    clusttaggen(i,:) =result{j,2}(1,i).add_clustTag;    
end
    
    clustTag{j,1}={clusttaggen};
end


clustName = cell(size(result, 1), 1);
for j = 1:size(result, 1)
    genSolution = result{j, 2};
    clustnamegen = zeros(numel(genSolution),num_clustName_columns);
for i = 1:numel(genSolution)
    clustnamegen(i,:) =result{j,2}(1,i).add_clustName;    
end
    clustName{j,1}={clustnamegen};
end

end

