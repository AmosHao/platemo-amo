function [pop, objVals, clustTag, clustName, centroid, Population] = ...
    EnvironmentalSelection(auxPop, pop, objVals, clustTag, clustName, ...
    centroid, selectedSize, objDim, varDim, Kmax, Problem)
% 环境选择函数（简化版，移除碰撞检测）
% 基于非支配排序和超体积指标进行选择
%
% 输入:
%   auxPop - 子代种群
%   pop - 父代种群
%   objVals - 父代目标值
%   clustTag - 聚类标签
%   clustName - 聚类名称
%   centroid - 聚类中心（目标空间）
%   selectedSize - 选择数量
%   objDim - 目标维度
%   varDim - 变量维度
%   Kmax - 最大聚类数
%   Problem - 问题对象
%
% 输出:
%   pop - 选择后的种群
%   objVals - 选择后的目标值
%   clustTag - 更新后的聚类标签
%   clustName - 更新后的聚类名称
%   centroid - 更新后的聚类中心
%   Population - 选择后的Population对象

    % 合并父代和子代
    if isempty(pop)
        popAll = auxPop;
    elseif isempty(auxPop)
        popAll = pop;
    else
        popAll = [pop; auxPop];
    end
    
    % 确保popAll不为空
    if isempty(popAll)
        % 如果都为空，返回空结果
        pop = [];
        objVals = [];
        clustTag = [];
        clustName = [];
        centroid = [];
        Population = [];
        return;
    end
    
    % 评估所有解
    allVals = Problem.Evaluation(popAll);
    objValsAll = allVals.objs;
    
    % 清理目标值，确保是实数
    objValsAll = real(objValsAll);
    objValsAll(isnan(objValsAll)) = 1e10;  % NaN替换为大数
    objValsAll(isinf(objValsAll)) = 1e10;  % Inf替换为大数
    
    % 调试：检查f1值的唯一性（每10代输出一次）
    persistent debugGen;
    if isempty(debugGen)
        debugGen = 0;
    end
    debugGen = debugGen + 1;
    if mod(debugGen, 10) == 1 && size(objValsAll, 1) > 0
        % 排除无效解（1e10）后统计
        validIdx = objValsAll(:, 1) < 1e9;  % 排除无效解
        validF1 = objValsAll(validIdx, 1);
        invalidCount = sum(~validIdx);
        
        if ~isempty(validF1)
            uniqueF1 = unique(round(validF1, 6));  % 保留6位小数
            fprintf('调试 [代 %d]: 总解数=%d, 有效解=%d, 无效解=%d, 唯一f1值数=%d, f1范围=[%.2f, %.2f], f1标准差=%.6f\n', ...
                debugGen, size(objValsAll, 1), length(validF1), invalidCount, length(uniqueF1), ...
                min(validF1), max(validF1), std(validF1));
            if length(uniqueF1) < length(validF1) * 0.5
                fprintf('警告: f1值多样性不足！只有%.1f%%的有效解有唯一的f1值\n', ...
                    100 * length(uniqueF1) / length(validF1));
            end
            if invalidCount > size(objValsAll, 1) * 0.1
                fprintf('警告: 无效解过多！%.1f%%的解是无效的（目标值=1e10）\n', ...
                    100 * invalidCount / size(objValsAll, 1));
            end
        else
            fprintf('调试 [代 %d]: 所有解都无效！\n', debugGen);
        end
    end
    
    % 更新聚类标签（新解标记为inf）
    tmpSize = size(auxPop, 1);
    clustTagAll = [clustTag; inf(tmpSize, 1)];
    
    % 非支配排序
    [FrontNo, ~] = NDSort(objValsAll, inf);
    
    % 构建参考点（改进：使用更合理的参考点，避免过度惩罚）
    % 使用非支配前沿的最大值，而不是所有解的最大值
    front1Indices = find(FrontNo == 1);
    if ~isempty(front1Indices)
        front1Max = max(objValsAll(front1Indices, :), [], 1);
        % 参考点：使用第一层的最大值，但不要太小，避免左下角解一直被保留
        refPoint = 1.15 * front1Max;  % 稍微扩大，确保有足够的淘汰空间
    else
        refPoint = max(objValsAll, [], 1);
        refPoint = 1.2 * refPoint;
    end
    
    % 按非支配层选择
    auxVals = [];
    auxPop = [];
    auxTag = [];
    auxSize = 0;
    auxAll = [];
    
    i = 1;
    FiSize = sum(FrontNo == i);
    maxFront = max(FrontNo(FrontNo ~= inf));
    if isempty(maxFront)
        maxFront = 1;
    end
    
    while auxSize < selectedSize && i <= maxFront
        frontIndices = find(FrontNo == i);
        if ~isempty(frontIndices)
            auxVals = [auxVals; objValsAll(frontIndices, 1:objDim)];
            auxPop = [auxPop; popAll(frontIndices, 1:varDim)];
            auxTag = [auxTag; clustTagAll(frontIndices, 1)];
            auxSize = auxSize + FiSize;
            auxAll = [auxAll, allVals(frontIndices)];
        end
        i = i + 1;
        if i <= maxFront
            FiSize = sum(FrontNo == i);
        else
            FiSize = 0;
        end
    end
    
    if i > 1
        FiSize = sum(FrontNo == i-1);
    else
        FiSize = 0;
    end
    
    % 如果超过选择数量，使用超体积指标删除多余解
    if i-1 ~= 1 && auxSize > selectedSize
        pop = auxPop;
        objVals = auxVals;
        clustTag = auxTag;
        PSize = auxSize;
        
        while PSize > selectedSize && PSize > 0
            if objDim == 2 && PSize >= 2
                % 2目标：使用拥挤距离（改进：处理f1值相同的情况）
                frontObjvs = objVals;
                fitV = zeros(PSize, 1);
                [frontObjvs, IX] = sortrows(frontObjvs, [1, 2]);  % 先按f1排序，f1相同时按f2排序
                
                % 检查f1值的唯一性
                uniqueF1 = unique(frontObjvs(:, 1));
                hasDuplicateF1 = length(uniqueF1) < PSize;
                
                if hasDuplicateF1
                    % 如果f1值有重复，优先删除f2值最接近其他解的解（减少多样性）
                    % 或者优先保留f2值差异最大的解
                    for i = 1:PSize
                        % 计算该解与其他解的f2距离
                        f2Distances = abs(frontObjvs(i, 2) - frontObjvs(:, 2));
                        f2Distances(i) = inf;  % 排除自己
                        minF2Dist = min(f2Distances);
                        
                        % 如果f1值相同，使用f2距离作为拥挤距离
                        sameF1Idx = find(abs(frontObjvs(:, 1) - frontObjvs(i, 1)) < 1e-10);
                        if length(sameF1Idx) > 1
                            % f1值相同，使用f2距离
                            fitV(i) = minF2Dist;
                        else
                            % f1值不同，使用标准拥挤距离
                            if i == 1
                                fitV(i) = (frontObjvs(2, 1) - frontObjvs(1, 1)) * ...
                                    (refPoint(1, 2) - frontObjvs(1, 2));
                            elseif i == PSize
                                fitV(i) = (refPoint(1, 1) - frontObjvs(PSize, 1)) * ...
                                    (frontObjvs(PSize-1, 2) - frontObjvs(PSize, 2));
                            else
                                fitV(i) = (frontObjvs(i+1, 1) - frontObjvs(i-1, 1)) * ...
                                    (frontObjvs(i-1, 2) - frontObjvs(i+1, 2));
                            end
                        end
                    end
                else
                    % 标准拥挤距离计算（f1值都不同）
                    if PSize == 2
                        fitV(IX(1)) = (frontObjvs(2, 1) - frontObjvs(1, 1)) * ...
                            (refPoint(1, 2) - frontObjvs(1, 2));
                        fitV(IX(2)) = (refPoint(1, 1) - frontObjvs(2, 1)) * ...
                            (frontObjvs(1, 2) - frontObjvs(2, 2));
                    else
                        fitV(IX(1)) = (frontObjvs(2, 1) - frontObjvs(1, 1)) * ...
                            (refPoint(1, 2) - frontObjvs(1, 2));
                        fitV(IX(2:PSize-1)) = (frontObjvs(3:PSize, 1) - frontObjvs(1:PSize-2, 1)) .* ...
                            abs(frontObjvs(1:PSize-2, 2) - frontObjvs(3:PSize, 2));
                        fitV(IX(PSize)) = (refPoint(1, 1) - frontObjvs(PSize, 1)) * ...
                            (frontObjvs(PSize-1, 2) - frontObjvs(PSize, 2));
                    end
                end
                
                % 如果拥挤距离为0或太小，使用f2距离作为备选
                fitV(fitV < 1e-10) = 0;
                if all(fitV == 0)
                    % 所有拥挤距离都为0，使用f2距离
                    for i = 1:PSize
                        f2Distances = abs(frontObjvs(i, 2) - frontObjvs(:, 2));
                        f2Distances(i) = inf;
                        fitV(i) = min(f2Distances);
                    end
                end
                
                % 改进：当f1值相同时，优先保留f2值差异大的解
                % 如果所有拥挤距离都很小，添加小的随机扰动，避免总是删除同一个解
                if max(fitV) < 1e-6
                    fitV = fitV + 1e-8 * rand(PSize, 1);  % 添加微小随机扰动
                end
                
                [~, loc] = min(fitV);
            else
                % 多目标：使用超体积贡献
                if PSize > 1
                    totalHV = CalHV(objVals, refPoint);
                    fitV = zeros(PSize, 1);
                    for j = 1:PSize
                        tmpObjVals = objVals;
                        tmpObjVals(j, :) = [];
                        tmpHV = CalHV(tmpObjVals, refPoint);
                        fitV(j) = totalHV - tmpHV;
                    end
                    [~, loc] = min(fitV);
                else
                    loc = 1;  % 如果只有一个解，删除它（虽然不应该发生）
                end
            end
            
            if loc > 0 && loc <= PSize
                objVals(loc, :) = [];
                pop(loc, :) = [];
                clustTag(loc, :) = [];
                if length(auxAll) >= loc
                    auxAll(loc) = [];
                end
                PSize = PSize - 1;
            else
                break;  % 如果索引无效，退出循环
            end
        end
    else
        pop = auxPop;
        objVals = auxVals;
        clustTag = auxTag;
    end
    
    % 决策空间多样性维护：去除重复个体，用批量随机新解替换（一次初始化，保证速度）
    [~, ia, ~] = unique(pop, 'rows');
    duplicateIndices = setdiff(1:size(pop, 1), ia);
    if ~isempty(duplicateIndices)
        numReplace = length(duplicateIndices);
        newPop = Problem.Initialization(numReplace);
        pop(duplicateIndices, :) = newPop.decs;
        objVals(duplicateIndices, :) = newPop.objs;
        clustTag(duplicateIndices) = inf;
        % 用 newPop 中已有解直接填 auxAll，避免逐次 Evaluation 拖慢迭代
        if numReplace == 1
            auxAll(duplicateIndices(1)) = newPop(1);
        else
            for j = 1:numReplace
                auxAll(duplicateIndices(j)) = newPop(j);
            end
        end
    end
    
    % 更新聚类信息（目标空间聚类）
    % 先为未分类的解分配聚类
    unclassified = find(clustTag == inf);
    for j = 1:length(unclassified)
        idx = unclassified(j);
        if idx > size(objVals, 1) || idx < 1
            continue;  % 跳过无效索引
        end
        
        % 获取当前解的目标值，确保是实数
        currentObj = objVals(idx, :);
        currentObj = real(currentObj);  % 确保是实数
        currentObj(isnan(currentObj)) = 0;  % 将NaN替换为0
        currentObj(isinf(currentObj)) = 1e10;  % 将Inf替换为大数
        
        if isempty(clustName) || isempty(centroid)
            clustTag(idx) = 1;
            clustName = 1;
            centroid = currentObj;
        else
            % 清理centroid，确保是实数
            centroid = real(centroid);
            centroid(isnan(centroid)) = 0;
            centroid(isinf(centroid)) = 1e10;
            
            % 基于目标空间距离分配到最近聚类
            if size(centroid, 1) > 0 && size(centroid, 2) == size(currentObj, 2)
                distances = pdist2(currentObj, centroid);
                [~, nearestClust] = min(distances);
                if nearestClust <= length(clustName)
                    clustTag(idx) = clustName(nearestClust);
                else
                    clustTag(idx) = clustName(1);  % 如果索引越界，使用第一个聚类
                end
            else
                % 如果维度不匹配，创建新聚类
                newName = max(clustName) + 1;
                clustTag(idx) = newName;
                clustName = [clustName; newName];
                centroid = [centroid; currentObj];
            end
        end
    end
    
    % 更新聚类信息（合并/分割聚类）
    [clustTag, clustName, centroid] = ObjectiveSpaceClustering(...
        objVals, clustTag, clustName, centroid, Kmax);
    
    % 构建Population对象
    Population = auxAll;
end


function HV = CalHV(PopObj, refPoint)
% 计算超体积（2目标问题的精确计算）
% 使用参考点作为边界

    [N, M] = size(PopObj);
    if N == 0
        HV = 0;
        return;
    end
    
    % 确保目标值为有限值
    PopObj = real(PopObj);
    PopObj(isnan(PopObj)) = 1e10;
    PopObj(isinf(PopObj)) = 1e10;
    refPoint = real(refPoint);
    refPoint(isnan(refPoint)) = 1e10;
    refPoint(isinf(refPoint)) = 1e10;
    
    % 对于2目标问题，使用精确计算
    if M == 2
        % 按第一个目标排序（升序）
        [sortedObj, ~] = sortrows(PopObj, 1);
        
        % 移除被支配的解（相对于参考点）
        validIdx = sortedObj(:, 1) < refPoint(1) & sortedObj(:, 2) < refPoint(2);
        if sum(validIdx) == 0
            HV = 0;
            return;
        end
        sortedObj = sortedObj(validIdx, :);
        N = size(sortedObj, 1);
        
        if N == 0
            HV = 0;
            return;
        end
        
        % 计算超体积（从参考点开始）
        HV = 0;
        prevX = refPoint(1);
        
        for i = 1:N
            % 当前解的x坐标
            currX = sortedObj(i, 1);
            % 宽度：从上一个x到当前x
            width = prevX - currX;
            
            % 高度：从当前y到参考点y（或到下一个更小的y）
            if i == N
                height = refPoint(2) - sortedObj(i, 2);
            else
                % 找到在当前x右侧且y更小的点
                rightPoints = sortedObj(i+1:N, :);
                if ~isempty(rightPoints)
                    minY = min(rightPoints(:, 2));
                    height = minY - sortedObj(i, 2);
                    if height < 0
                        height = refPoint(2) - sortedObj(i, 2);
                    end
                else
                    height = refPoint(2) - sortedObj(i, 2);
                end
            end
            
            if width > 0 && height > 0
                HV = HV + width * height;
            end
            
            prevX = currX;
        end
    else
        % 多目标：使用简化计算（每个解的贡献）
        HV = 0;
        for i = 1:N
            % 计算该解到参考点的超体积贡献
            diff = refPoint - PopObj(i, :);
            diff(diff < 0) = 0;  % 只考虑正的部分
            vol = prod(diff);
            HV = HV + vol;
        end
    end
end
