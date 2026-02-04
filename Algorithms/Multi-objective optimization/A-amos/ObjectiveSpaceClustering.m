function [clustTag, clustName, centroid] = ObjectiveSpaceClustering(objVals, clustTag, clustName, centroid, Kmax)
% 目标空间聚类函数
% 在目标空间对解进行聚类，用于指导邻域搜索
%
% 输入:
%   objVals - 目标值矩阵 (N x M)
%   clustTag - 当前聚类标签 (N x 1)
%   clustName - 聚类名称列表
%   centroid - 聚类中心（目标空间）(K x M)
%   Kmax - 最大聚类数
%
% 输出:
%   clustTag - 更新后的聚类标签
%   clustName - 更新后的聚类名称列表
%   centroid - 更新后的聚类中心（目标空间）

    N = size(objVals, 1);
    nClust = size(clustName, 1);
    
    % 清理数据，确保是实数
    objVals = real(objVals);
    objVals(isnan(objVals)) = 0;
    objVals(isinf(objVals)) = 1e10;
    
    if ~isempty(centroid)
        centroid = real(centroid);
        centroid(isnan(centroid)) = 0;
        centroid(isinf(centroid)) = 1e10;
    end
    
    % 如果聚类数超过Kmax，合并最近的聚类
    while nClust > Kmax && nClust > 1
        % 计算聚类中心之间的距离（在目标空间）
        if size(centroid, 1) < 2
            break;  % 如果聚类数少于2，无法合并
        end
        distances = pdist(centroid);
        disMatrix = squareform(distances);
        for i = 1:nClust
            disMatrix(i, i) = inf;
        end
        
        % 找到两个最近的聚类
        [minV, ic] = min(disMatrix(:));
        [ir, ic] = ind2sub(size(disMatrix), ic);
        
        % 合并聚类
        if isempty(clustName)
            newName = 1;
        else
            newName = max(clustName) + 1;
        end
        % 计算新聚类中心（目标空间均值）
        members1 = objVals(clustTag == clustName(ir), :);
        members2 = objVals(clustTag == clustName(ic), :);
        if ~isempty(members1) && ~isempty(members2)
            newCentroid = mean([members1; members2], 1);
        elseif ~isempty(members1)
            newCentroid = mean(members1, 1);
        elseif ~isempty(members2)
            newCentroid = mean(members2, 1);
        else
            newCentroid = centroid(ir, :);  % 如果都为空，使用第一个聚类中心
        end
        newCentroid = real(newCentroid);
        newCentroid(isnan(newCentroid)) = 0;
        newCentroid(isinf(newCentroid)) = 1e10;
        
        % 更新聚类标签
        clustTag(clustTag == clustName(ir)) = newName;
        clustTag(clustTag == clustName(ic)) = newName;
        
        % 更新聚类名称和中心
        centroid([ir; ic], :) = [];
        clustName([ir; ic], :) = [];
        clustName = [clustName; newName];
        centroid = [centroid; newCentroid];
        
        nClust = size(clustName, 1);
    end
    
    % 如果聚类数少于Kmax，为未分类的解创建新聚类
    unclassified = find(clustTag == inf | clustTag == 0);
    if ~isempty(unclassified) && ~isempty(clustName)
        for i = 1:length(unclassified)
            idx = unclassified(i);
            if idx > 0 && idx <= N
                newName = max(clustName) + 1;
                clustTag(idx) = newName;
                clustName = [clustName; newName];
                newCentroid = real(objVals(idx, :));
                newCentroid(isnan(newCentroid)) = 0;
                newCentroid(isinf(newCentroid)) = 1e10;
                centroid = [centroid; newCentroid];
            end
        end
    elseif ~isempty(unclassified) && isempty(clustName)
        % 如果没有任何聚类，创建第一个聚类
        idx = unclassified(1);
        if idx <= size(objVals, 1) && idx >= 1
            clustTag(idx) = 1;
            clustName = 1;
            newCentroid = real(objVals(idx, :));
            newCentroid(isnan(newCentroid)) = 0;
            newCentroid(isinf(newCentroid)) = 1e10;
            centroid = newCentroid;
        end
    end
end
