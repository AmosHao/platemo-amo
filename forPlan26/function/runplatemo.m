% [result,Obj]=runplatemo1(@huigui,@SA,100,1,3,100000,20);
% function [result,Obj]=runplatemo1(problem,algorithm,N,M,D,maxFE,s)


% [result,igd,hv]=runplatemo(@GLT1,@test1GOCEA_kmax5_hvtest,100,10000,100);
function [Obj,igd,hv,igd_end,hv_end]=runplatemo(problem,algorithm,N,M,D,maxFE,s)
platemo('problem',problem,'algorithm',algorithm,'N',N,'M',M,'D',D,'maxFE',maxFE,'save',s);
result = evalin('base', 'result');
% 获取最后一行最后一列的 SOLUTION 结构
lastSolution = result{end, end};
% 获取每个solution的obj值的列数
num_obj_columns = numel(lastSolution(1).obj);
% 提取每个1*1的solution的obj值
Obj = zeros(numel(lastSolution),num_obj_columns); % 初始化存储obj值的数组
for idx = 1:numel(lastSolution)
    solution = lastSolution(idx);
    obj = solution.obj; % 假设obj字段存储了obj的值
    Obj(idx,:) = obj;
end
% pro = problem('N',N,'M',M,'D',D,'maxFE',maxFE);
% igd_end=pro.CalMetric('IGD',result{end});
% hv_end=pro.CalMetric('HV',result{end});
% % 提取每代的IGD值
% igd = zeros(size(result, 1), 1);
% hv = zeros(size(result, 1), 1);
% for i = 1:size(result, 1)
%     igd(i) = pro.CalMetric('IGD', result{i, end});
%     hv(i) = pro.CalMetric('HV', result{i, end});
% end
end

