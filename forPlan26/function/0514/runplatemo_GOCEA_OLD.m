%[Obj,igd,hv,igd_end,hv_end]=runplatemo(@GLT1,@test1GOCEAmy,5,0.9,0.5,0.1,100,10000,51);
function [Obj,igd,hv,igd_end,hv_end]=runplatemo_GOCEA_OLD(problem,algorithm,K,BETA,F,CR,N,maxFE,s,vardim,M)
platemo('problem',problem,'algorithm',{algorithm,K,BETA,F,CR},'D',vardim,'N',N,'M',M,'maxFE',maxFE,'save',s);
result = evalin('base', 'result');

lastSolution = result{end, end};% 获取最后一行最后一列的 SOLUTION 结构

num_obj_columns = numel(lastSolution(1).obj);% 获取每个solution的obj值的列数

Obj = zeros(numel(lastSolution),num_obj_columns); % 初始化存储obj值的数组。
for idx = 1:numel(lastSolution)
    solution = lastSolution(idx);
    obj = solution.obj; % 提取每个1*1的solution的obj值
    Obj(idx,:) = obj;
end
pro = problem('N',N,'M',M,'D',vardim,'maxFE',maxFE);
igd_end=pro.CalMetric('IGD',result{end});
hv_end=pro.CalMetric('HV',result{end});
% 提取每代的IGD值
igd = zeros(size(result, 1), 1);
hv = zeros(size(result, 1), 1);
for i = 1:size(result, 1)
    igd(i) = pro.CalMetric('IGD', result{i, end});
    hv(i) = pro.CalMetric('HV', result{i, end});
end
end

