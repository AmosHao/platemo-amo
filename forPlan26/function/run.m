clear
problem={@LZ1,@LZ2};
[Obj,igd,hv,igd_end,hv_end]=runplatemo(problem{1},@test1GOCEA_kmax5_hvtest,100,10,2,60000,20);
% problem=GLT1;
% algorithm=@test1GOCEA_kmax5_hvtest;
% N=100;
% maxFE=10000;
% save=100;
function [Obj,igd,hv,igd_end,hv_end]=runplatemo(problem,algorithm,N,D,M,maxFE,save)
platemo('problem',problem,'algorithm',algorithm,'N',N,'D',D,'M',M,maxFE',maxFE,'save',save,'outputFcn');
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
% % 从 SOLUTION 结构中提取 obj
% obj = lastSolution.obj;
pro = problem();
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

% pro = problem();
% igd=pro.CalMetric('IGD',result{end});
% hv=pro.CalMetric('HV',result{end});

% [i, h] = runplatemo(@GLT1, @test1GOCEA_kmax5_hvtest, 100, 10000);
% function [igd,hv] = runplatemo(problem, algorithm, N, maxFE)
% % 调用 platemo 函数，传递问题函数、算法函数、N和maxFE
%     platemo_result = platemo('problem', problem, 'algorithm', algorithm, ...
%         'N', N, 'maxFE', maxFE, 'save', -10, 'outputFcn', []);
% 
%     % 获取 platemo 结果
%     last_result = platemo_result{end};
% 
%     % 通过问题函数获取问题对象
%     pro = problem_function();
% 
%     % 使用问题对象执行操作，计算IGD和HV
%     igd = pro.CalMetric('IGD', last_result);
%     hv = pro.CalMetric('HV', last_result);
% end