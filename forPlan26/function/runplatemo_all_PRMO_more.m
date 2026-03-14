%platemo('problem',@IMMOEA,'algorithm',@GLT1,'D',3,'N',100,'maxFE',10000,'save',20);
clear
% dotss=[500,500,800;
%     7000,3000,800;
%     9500,6000,800;
%     6050,9200,800;
%     1800,8950,800;
%     2800,4050,800;
%     6000,6100,800];
% 客户点坐标矩阵
        dotss = [
            5500, 5000, 800;   % 配送中心（第1行）
             500,  500, 800;   % 客户点1
            3800, 2800, 800;   % 客户点2
            7000, 3000, 800;   % 客户点3
            9000, 4000, 800;   % 客户点4
            2800, 4500, 800;   % 客户点5
            9500, 6000, 800;   % 客户点6
             500, 7000, 800;   % 客户点7
            6000, 6100, 800;   % 客户点8
            1800, 8950, 800;   % 客户点9
            6050, 9200, 800    % 客户点10
        ];

combines=nchoosek(1:size(dotss,1),2);
combinesnum=size(combines,1);
for k =1:combinesnum
    % if any(combines(k,:)==7)
    % dots=[dotss(combines(k,1),:),dotss(combines(k,2),:),];
    % dot1=combines(k,1);
    % dot2=combines(k,2);
    if k>4
    dots=[dotss(combines(k,1),:),dotss(combines(k,2),:),];
    dot1=combines(k,1);
    dot2=combines(k,2);
% [validtrait_end,validtrait,clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(@planner_travel_maxhjd_newobj,@GOCEAfortravel_singlepath_5,100,60000,20,3,dots,dot1,dot2);
[validtrait_end,validtrait,Dec,Obj,Dec_end,Obj_end,hv,hv_end]=runplatemo_all1(@planner_travel_maxhjd_newobj,@PRMO,100,60000,20,3,dots,dot1,dot2);
% [clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(@planner_travel_maxhjd02,@RMDE,100,10000,20,3,dots,dot1,dot2);
    end
end
% function [validtrait_end,validtrait,clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(problem,algorithm,N,maxFE,s,M,dots,dot1,dot2)
% function [validtrait_end,validtrait,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all1(problem,algorithm,N,maxFE,s,M,dots,dot1,dot2)

% function [validtrait_end,validtrait,Dec,Obj,Dec_end,Obj_end,hv,hv_end,pd,pd_end]=runplatemo_all_singlepath_PRMO(problem,algorithm,N,maxFE,s,M,dots)
% function [Dec,Obj,Dec_end,Obj_end,hv,hv_end,pd,pd_end]=runplatemo_all_singlepath_PRMO(problem,algorithm,N,maxFE,s,M,dots)

function [validtrait_end,validtrait,Dec,Obj,Dec_end,Obj_end,hv,hv_end]=runplatemo_all1(problem,algorithm,N,maxFE,s,M,dots,dot1,dot2)

platemo('problem',problem,'algorithm',algorithm,'N',N,'M',M,'maxFE',maxFE,'save',s,'dots',dots);
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
validtrait_end = zeros(numel(lastSolution),1);
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
    validtrait = solution.validtrait;
    validtrait_end(idx,:) = validtrait;
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

pro = problem('M',M,'dots',dots);
% igd_end=pro.CalMetric('IGD',result{end});
hv_end=pro.CalMetric('HV',result{end});
pd_end=pro.CalMetric('PD',result{end});
% 提取每代的IGD值
% igd = zeros(size(result, 1), 1);
hv = zeros(size(result, 1), 1);
pd = zeros(size(result, 1), 1);
for i = 1:size(result, 1)
    % igd(i) = pro.CalMetric('IGD', result{i, end});
    hv(i) = pro.CalMetric('HV', result{i, end});
    pd(i) = pro.CalMetric('PD', result{i, end});
end

objgen = zeros(numel(lastSolution), num_obj_columns);
Obj = cell(size(result, 1), 1);
for j = 1:size(result, 1)
for i = 1:numel(result{j, end})
    % objgen = zeros(numel(result{j, end}), num_obj_columns);
    objgen(i,:) =result{j,2}(1,i).obj;    
end
    
    Obj{j,1}={objgen};
end



popgen = zeros(numel(lastSolution), num_dec_columns);
Dec = cell(size(result, 1), 1);
for j = 1:size(result, 1)
for i = 1:numel(result{j, end})
    % popgen = zeros(numel(result{j, end}), num_dec_columns);
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
% validtraitgen = zeros(numel(result{j, end}), 1);
validtrait = cell(size(result, 1), 1);
for j = 1:size(result, 1)
    validtraitgen = zeros(numel(result{j, end}), 1);
for i = 1:numel(result{j, end})
    validtraitgen(i,:) =result{j,2}(1,i).validtrait;    
end

    validtrait{j,1}={validtraitgen};
end


% end

% 指定要保存的Excel文件名
% filename = 'maxhjddot2_disselect_initialselect_objduiying01.xlsx';
% 将 dot1 和 dot2 转换为字符串，用于文件名
dot1_str = num2str(dot1);
dot2_str = num2str(dot2);

folder = '/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/PlatEMO-master-using/order/Data0825';          % ← 换成你的目标文件夹

% 创建唯一的文件名
filename = fullfile(folder, sprintf('0825_dot1_%s_dot2_%s.xlsx', dot1_str, dot2_str));

% 创建唯一的文件名
% filename = sprintf('0825_dot1_%s_dot2_%s.xlsx', dot1_str, dot2_str);

% 将 obj 矩阵保存到 'Sheet1'
writematrix(Obj_end, filename, 'Sheet', 'obj_data');

% 将 dec 矩阵保存到 'Sheet2'
writematrix(Dec_end, filename, 'Sheet', 'dec_data');

writematrix(validtrait_end, filename, 'Sheet', 'validtrait_data');

% writematrix(bj_end, filename, 'Sheet', 'bj_data');

writematrix(hv, filename, 'Sheet', 'hv_data');
end


