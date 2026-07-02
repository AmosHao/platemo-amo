% 与 runplatemo_all_singlepath_PRMO_class 相同，用于 fangzhen 问题与 PRMO_fangzhen
% 搭配：@planner_travel_maxhjd_newobj_fangzhen, @PRMO_fangzhen
function [validtrait_end,validtrait,clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,hv,hv_end,pd,pd_end]=runplatemo_all_singlepath_PRMO_fangzhen(problem,algorithm,N,maxFE,s,M,dots)
% 注意：不能在 parfor worker 中使用 evalin('base',...)（会触发 Transparency violation）。
% 这里直接构造 Problem，并用 Algorithm.Solve 获取 Algorithm.result（与原 platemo 保存逻辑一致）。
[Problem, Algorithm] = platemo_build_problem_alg(problem, algorithm, N, maxFE, M, dots, s);
Algorithm.Solve(Problem);
result = Algorithm.result;
lastSolution = result{end, end};
a=numel(lastSolution);
num_obj_columns = numel(lastSolution(1).obj);
num_dec_columns = numel(lastSolution(1).dec);
num_clustName_columns = numel(lastSolution(1).add_clustName);
num_clustTag_columns = numel(lastSolution(1).add_clustTag);

Obj_end = zeros(numel(lastSolution),num_obj_columns);
Dec_end = zeros(numel(lastSolution),num_dec_columns);
clustTag_end = zeros(numel(lastSolution),num_clustTag_columns);
clustName_end = zeros(numel(lastSolution),num_clustName_columns);
validtrait_end = zeros(numel(lastSolution),1);
for idx = 1:numel(lastSolution)
    solution = lastSolution(idx);
    Obj = solution.obj;
    Obj_end(idx,:) = Obj;
    dec = solution.dec;
    Dec_end(idx,:) = dec;
    add_clustTag = solution.add_clustTag;
    clustTag_end(idx,:) = add_clustTag;
    add_clustName = solution.add_clustName;
    clustName_end(idx,:) = add_clustName;
    validtrait = solution.validtrait;
    if isempty(validtrait)
        validtrait_end(idx,:) = 0;
    else
        validtrait_end(idx,:) = validtrait;
    end
end

pro = problem('M',M,'dots',dots);
% parpool worker（包含 parfor/spmd lab）上计算某些指标可能会触发 Transparency violation；指标不影响 edge.mat 的核心 Dec/Obj
skipMetrics = local_should_skip_metrics_for_parallel();

hv_end = NaN; pd_end = NaN;
hv = NaN(size(result, 1), 1);
pd = NaN(size(result, 1), 1);

if ~skipMetrics
try
    hv_end = pro.CalMetric('HV',result{end});
catch
end
try
    pd_end = pro.CalMetric('PD',result{end});
catch
end
for i = 1:size(result, 1)
    try
        hv(i) = pro.CalMetric('HV', result{i, end});
    catch
        hv(i) = NaN;
    end
    try
        pd(i) = pro.CalMetric('PD', result{i, end});
    catch
        pd(i) = NaN;
    end
end
end

objgen = zeros(numel(lastSolution), num_obj_columns);
Obj = cell(size(result, 1), 1);
for j = 1:size(result, 1)
for i = 1:numel(result{j, end})
    objgen(i,:) =result{j,2}(1,i).obj;    
end
    Obj{j,1}={objgen};
end

popgen = zeros(numel(lastSolution), num_dec_columns);
Dec = cell(size(result, 1), 1);
for j = 1:size(result, 1)
for i = 1:numel(result{j, end})
    popgen(i,:) =result{j,2}(1,i).dec;    
end
    Dec{j,1}={popgen};
end

clusttaggen = zeros(numel(lastSolution), num_clustTag_columns);
clustTag = cell(size(result, 1), 1);
for j = 1:size(result, 1)
for i = 1:numel(result{j, end})
    tmpTag = result{j,2}(1,i).add_clustTag;
    if isempty(tmpTag)
        % 某些代/个体可能尚未写入该字段（或写入空），用 0 填充避免维度错误
        clusttaggen(i,:) = 0;
    else
        clusttaggen(i,:) = tmpTag;
    end
end
    clustTag{j,1}={clusttaggen};
end

clustnamegen = zeros(numel(lastSolution),num_clustName_columns);
clustName = cell(size(result, 1), 1);
for j = 1:size(result, 1)
for i = 1:numel(result{j, end})
    tmpName = result{j,2}(1,i).add_clustName;
    if isempty(tmpName)
        clustnamegen(i,:) = 0;
    else
        clustnamegen(i,:) = tmpName;
    end
end
    clustName{j,1}={clustnamegen};
end

validtrait = cell(size(result, 1), 1);
for j = 1:size(result, 1)
    validtraitgen = zeros(numel(result{j, end}), 1);
    for i = 1:numel(result{j, end})
        tmpValid = result{j,2}(1,i).validtrait;
        if isempty(tmpValid)
            validtraitgen(i,:) = 0;
        else
            validtraitgen(i,:) = tmpValid;
        end
    end
    validtrait{j,1}={validtraitgen};
end

end

function tf = local_should_skip_metrics_for_parallel()
% MATLAB 并行环境里：getCurrentTask 在 spmd worker 往往为空；用 labindex 更可靠
    tf = false;
    try
        if ~isempty(getCurrentTask())
            tf = true;
            return;
        end
    catch
    end
    try
        % labindex 仅在 spmd 内可用；外部调用会报错，捕获即可
        li = labindex;
        if ~isempty(li)
            tf = true;
        end
    catch
    end
end

function [Problem, Algorithm] = platemo_build_problem_alg(problem, algorithm, N, maxFE, M, dots, saveCount)
% 复制 platemo.m 的核心装配逻辑（去掉对 base workspace 的依赖）。
    if verLessThan('matlab','9.4')
        error('Fail to use PlatEMO since the version for MATLAB is lower than R2018a.');
    end
    Setting = {'problem',problem,'algorithm',algorithm,'N',N,'M',M,'maxFE',maxFE,'dots',dots,'save',saveCount};
    [PRO,input] = platemo_getSetting(Setting);
    Problem     = PRO(input{:});
    [ALG,input] = platemo_getSetting(Setting,Problem);
    Algorithm   = ALG(input{:});
end

function [name,Setting] = platemo_getSetting(Setting,Pro)
% 复制自 PlatEMO/platemo.m 的 getSetting（局部副本）。
    isStr = find(cellfun(@ischar,Setting(1:end-1))&~cellfun(@isempty,Setting(2:end)));
    if nargin > 1
        index = isStr(find(strcmp(Setting(isStr),'algorithm'),1)) + 1;
        if isempty(index)
            names = {@BSPGA,@GA,@SACOSO,@GA;@PMMOEA,@NSGAIII,@KRVEA,@NSGAIII;@RVEA,@RVEA,@CSEA,@RVEA};
            name  = names{find([Pro.M<2,Pro.M<4,1],1),find([all(Pro.encoding==4),any(Pro.encoding>2),Pro.maxFE<=1000&Pro.D<=10,1],1)};
        elseif iscell(Setting{index})
            name    = Setting{index}{1};
            Setting = [Setting,{'parameter'},{Setting{index}(2:end)}];
        else
            name = Setting{index};
        end
    else
        index = isStr(find(strcmp(Setting(isStr),'problem'),1)) + 1;
        if isempty(index)
            name = @UserProblem;
        elseif iscell(Setting{index})
            name    = Setting{index}{1};
            Setting = [Setting,{'parameter'},{Setting{index}(2:end)}];
        else
            name = Setting{index};
        end
    end
end
