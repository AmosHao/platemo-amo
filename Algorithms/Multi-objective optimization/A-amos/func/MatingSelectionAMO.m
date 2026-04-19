function [MatingPool1,MatingPool2,state] = MatingSelectionAMO(Population,N,varargin)
% MatingSelectionAMO - 目标空间聚类引导的交配选择（可嵌入任意 MOEA）
%
% 用法：
%   [p1,p2] = MatingSelectionAMO(Population,N)
%   [p1,p2,state] = MatingSelectionAMO(Population,N,'Kmax',20,'state',state,'p1',p1)
%
% 输入：
%   Population : INDIVIDUAL 数组（需有 .objs；若需要沿用某算法的父代1选择，可传入 p1）
%   N          : 需要产生的配对数量（通常等于 Problem.N）
%
% 参数（键值对，可选）：
%   Kmax  : 最大聚类数（默认 20）
%   state : 结构体，允许跨代复用 clustTag/clustName/centroid（可为空）
%   p1    : 已生成的父代1索引向量（长度 N）。若不提供，则默认随机从 1..|Pop| 采样
%
% 输出：
%   MatingPool1 : 父代1索引（1×N）
%   MatingPool2 : 父代2索引（1×N），优先跨簇选择；无跨簇候选时退化为全局随机（除自身）
%   state       : 更新后的聚类状态（clustTag/clustName/centroid）
%
% 说明：
% - 本函数只负责“聚类交配选择”。你可以：
%   1) 保持原算法的父代1选择（锦标赛/邻域等），把结果传入 'p1'
%   2) 或者不传 p1，让它随机产生 p1
% - 本实现采用“完全增量聚类”：
%   * 若 state.centroid 不存在/维度不匹配，则初始化一次聚类得到 centroid
%   * 后续每代：用上代 centroid 进行一次“分配→更新”，不再从头重聚类

    % ---- defaults ----
    Kmax = 20;
    state = struct();
    p1Provided = false;
    MatingPool1 = [];

    % ---- parse args ----
    if ~isempty(varargin)
        for i = 1:2:length(varargin)
            key = lower(string(varargin{i}));
            val = varargin{i+1};
            switch key
                case "kmax"
                    Kmax = val;
                case "state"
                    state = val;
                case {"p1","matingpool1"}
                    MatingPool1 = val;
                    p1Provided = true;
            end
        end
    end

    popSize = length(Population);
    if popSize < 1
        MatingPool1 = [];
        MatingPool2 = [];
        return;
    end

    if ~p1Provided || isempty(MatingPool1)
        MatingPool1 = randi(popSize,1,N);
    else
        MatingPool1 = reshape(MatingPool1,1,[]);
        if length(MatingPool1) ~= N
            MatingPool1 = MatingPool1(1:min(N,end));
            if length(MatingPool1) < N
                MatingPool1 = [MatingPool1, randi(popSize,1,N-length(MatingPool1))];
            end
        end
    end

    % ---- incremental clustering in objective space ----
    K = max(1,min(Kmax,popSize));
    Obj = real(Population.objs);
    Obj(isnan(Obj)) = 0;
    Obj(isinf(Obj)) = 1e10;
    [~,M] = size(Obj);

    needInit = ~isfield(state,'centroid') || isempty(state.centroid) || size(state.centroid,1)~=K || size(state.centroid,2)~=M;
    if needInit
        % 初始化一次：用简单 kmeans 得到初始中心与簇标签
        [state.clustTag,state.centroid] = KMeansObjectiveClustering(Obj,K);
        state.clustName = (1:K)';
    else
        % 分配：按上代中心给当前种群分簇
        centroid = state.centroid;
        d = zeros(popSize,K);
        for j = 1:K
            diff = Obj - centroid(j,:);
            d(:,j) = sum(diff.^2,2);
        end
        [minDist,tag] = min(d,[],2);

        % 更新：每个簇中心取当前簇内均值；空簇用“最远点”重置
        newCentroid = centroid;
        for j = 1:K
            members = (tag == j);
            if any(members)
                newCentroid(j,:) = mean(Obj(members,:),1);
            else
                % 选择离其所属中心最远的点来填补空簇（提高多样性）
                [~,idxFar] = max(minDist);
                newCentroid(j,:) = Obj(idxFar,:);
                tag(idxFar) = j;
                minDist(idxFar) = -inf;
            end
        end
        state.clustTag = tag;
        state.centroid = newCentroid;
        state.clustName = (1:K)';
    end

    % ---- cross-cluster mate selection ----
    MatingPool2 = zeros(1,N);
    for k = 1:N
        idx1 = MatingPool1(k);
        tag1 = state.clustTag(idx1);
        otherIdx = find(state.clustTag ~= tag1);
        otherIdx = setdiff(otherIdx, idx1);
        if isempty(otherIdx)
            otherIdx = setdiff(1:popSize, idx1);
        end
        if isempty(otherIdx)
            idx2 = idx1;
        else
            idx2 = otherIdx(randi(length(otherIdx)));
        end
        MatingPool2(k) = idx2;
    end
end

function [clustTag,centroid] = KMeansObjectiveClustering(objVals,K)
% Simple k-means clustering in objective space (no toolbox dependency).
% objVals: N×M, K: number of clusters (1..N)

    objVals = real(objVals);
    objVals(isnan(objVals)) = 0;
    objVals(isinf(objVals)) = 1e10;

    [N,M] = size(objVals);
    K = max(1,min(K,N));

    % Initialize centroids by random distinct samples
    if K == 1
        clustTag = ones(N,1);
        centroid = mean(objVals,1);
        return;
    end
    initIdx  = randperm(N,K);
    centroid = objVals(initIdx,:);

    clustTag = ones(N,1);
    maxIter  = 10;
    for iter = 1:maxIter
        % Assign step: nearest centroid (squared Euclidean)
        d = zeros(N,K);
        for j = 1:K
            diff = objVals - centroid(j,:);
            d(:,j) = sum(diff.^2,2);
        end
        [~,newTag] = min(d,[],2);

        if iter > 1 && all(newTag == clustTag)
            break;
        end
        clustTag = newTag;

        % Update step: recompute centroid; re-seed empty clusters
        for j = 1:K
            members = (clustTag == j);
            if any(members)
                centroid(j,:) = mean(objVals(members,:),1);
            else
                centroid(j,:) = objVals(randi(N),:);
            end
        end
    end
end

