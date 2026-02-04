function offspring = DiscreteMutation(individual, m)
% 离散序列变异操作
% 支持多种变异策略：交换、插入、段内重排、段间移动
%
% 输入:
%   individual - 待变异的个体
%   m - 无人机数量（段数）
%
% 输出:
%   offspring - 变异后的个体

    mutationType = randi(4);  % 随机选择变异类型
    
    switch mutationType
        case 1
            % 交换变异：随机选择两个非0位置交换
            offspring = SwapMutation(individual);
        case 2
            % 插入变异：随机选择一个点，插入到新位置
            offspring = InsertMutation(individual);
        case 3
            % 段内重排：在保持0分隔符的前提下，重排某段内的点
            offspring = SegmentReorderMutation(individual, m);
        case 4
            % 段间移动：将某个点从一个段移动到另一个段
            offspring = SegmentMoveMutation(individual, m);
        otherwise
            offspring = individual;
    end
end

function offspring = SwapMutation(individual)
% 交换变异
    nonZeroIdx = find(individual ~= 0);
    if length(nonZeroIdx) < 2
        offspring = individual;
        return;
    end
    
    idx = randperm(length(nonZeroIdx), 2);
    offspring = individual;
    offspring(nonZeroIdx(idx(1))) = individual(nonZeroIdx(idx(2)));
    offspring(nonZeroIdx(idx(2))) = individual(nonZeroIdx(idx(1)));
end

function offspring = InsertMutation(individual)
% 插入变异
    originalDim = length(individual);
    nonZeroIdx = find(individual ~= 0);
    if length(nonZeroIdx) < 2
        offspring = individual;
        return;
    end
    
    % 随机选择一个点
    fromIdx = nonZeroIdx(randi(length(nonZeroIdx)));
    val = individual(fromIdx);
    
    % 随机选择插入位置（非0位置，但不能是同一个位置）
    validToIdx = nonZeroIdx(nonZeroIdx ~= fromIdx);
    if isempty(validToIdx)
        offspring = individual;
        return;
    end
    toIdx = validToIdx(randi(length(validToIdx)));
    
    offspring = individual;
    offspring(fromIdx) = [];  % 删除元素后，数组长度变成 originalDim - 1
    
    % 找到新的插入位置（考虑删除后的数组）
    % 删除 fromIdx 后，所有在 fromIdx 之后的索引都要减1
    if toIdx > fromIdx
        newToIdx = toIdx - 1;  % 因为删除了 fromIdx，所以 toIdx 需要减1
    else
        newToIdx = toIdx;  % toIdx 在 fromIdx 之前，不受影响
    end
    
    % 确保 newToIdx 在有效范围内 [1, originalDim-1]
    newToIdx = max(1, min(newToIdx, originalDim - 1));
    
    % 插入：确保索引不超出范围
    if newToIdx >= length(offspring)
        % 如果插入位置在末尾或之后，直接追加
        offspring = [offspring, val];
    elseif newToIdx <= 0
        % 如果插入位置在开头或之前，在开头插入
        offspring = [val, offspring];
    else
        % 正常插入
        offspring = [offspring(1:newToIdx), val, offspring(newToIdx+1:end)];
    end
    
    % 确保维度一致
    if length(offspring) ~= originalDim
        if length(offspring) < originalDim
            offspring = [offspring, zeros(1, originalDim - length(offspring))];
        else
            offspring = offspring(1:originalDim);
        end
    end
end

function offspring = SegmentReorderMutation(individual, m)
% 段内重排变异
    originalDim = length(individual);
    idx0 = find(individual == 0);
    
    % 确保有 m-1 个0分隔符
    if length(idx0) < m - 1
        % 如果0的数量不足，在末尾补充0
        numZerosNeeded = (m - 1) - length(idx0);
        individual = [individual, zeros(1, numZerosNeeded)];
        idx0 = find(individual == 0);
    elseif length(idx0) > m - 1
        % 如果0的数量过多，只保留前m-1个
        idx0 = idx0(1:m-1);
    end
    
    idx0 = [0 idx0 length(individual)+1];
    
    % 随机选择一个段（确保不越界）
    segIdx = randi(min(m, length(idx0)-1));
    if segIdx < length(idx0)
        segRange = (idx0(segIdx)+1):(idx0(segIdx+1)-1);
    else
        segRange = (idx0(segIdx)+1):(idx0(end)-1);
    end
    
    if length(segRange) < 2 || isempty(segRange)
        offspring = individual;
    else
        % 重排该段内的点
        offspring = individual;
        segVals = offspring(segRange);
        segVals = segVals(randperm(length(segVals)));
        offspring(segRange) = segVals;
    end
    
    % 确保维度一致
    if length(offspring) ~= originalDim
        if length(offspring) < originalDim
            offspring = [offspring, zeros(1, originalDim - length(offspring))];
        else
            offspring = offspring(1:originalDim);
        end
    end
end

function offspring = SegmentMoveMutation(individual, m)
% 段间移动变异
    originalDim = length(individual);
    idx0 = find(individual == 0);
    
    % 确保有 m-1 个0分隔符
    if length(idx0) < m - 1
        % 如果0的数量不足，在末尾补充0
        numZerosNeeded = (m - 1) - length(idx0);
        individual = [individual, zeros(1, numZerosNeeded)];
        idx0 = find(individual == 0);
    elseif length(idx0) > m - 1
        % 如果0的数量过多，只保留前m-1个
        idx0 = idx0(1:m-1);
    end
    
    idx0 = [0 idx0 length(individual)+1];
    
    % 随机选择源段和目标段（确保不越界）
    maxSeg = min(m, length(idx0)-1);
    fromSeg = randi(maxSeg);
    toSeg = randi(maxSeg);
    
    if fromSeg == toSeg
        offspring = individual;
        % 确保维度一致
        if length(offspring) ~= originalDim
            if length(offspring) < originalDim
                offspring = [offspring, zeros(1, originalDim - length(offspring))];
            else
                offspring = offspring(1:originalDim);
            end
        end
        return;
    end
    
    % 安全访问段范围
    if fromSeg < length(idx0)
        if fromSeg+1 <= length(idx0)
            fromRange = (idx0(fromSeg)+1):(idx0(fromSeg+1)-1);
        else
            fromRange = (idx0(fromSeg)+1):(idx0(end)-1);
        end
    else
        fromRange = [];
    end
    
    if toSeg < length(idx0)
        if toSeg+1 <= length(idx0)
            toRange = (idx0(toSeg)+1):(idx0(toSeg+1)-1);
        else
            toRange = (idx0(toSeg)+1):(idx0(end)-1);
        end
    else
        toRange = [];
    end
    
    if isempty(fromRange)
        offspring = individual;
        % 确保维度一致
        if length(offspring) ~= originalDim
            if length(offspring) < originalDim
                offspring = [offspring, zeros(1, originalDim - length(offspring))];
            else
                offspring = offspring(1:originalDim);
            end
        end
        return;
    end
    
    % 从源段随机选择一个点
    pointIdx = fromRange(randi(length(fromRange)));
    pointVal = individual(pointIdx);
    
    % 移动到目标段
    offspring = individual;
    offspring(pointIdx) = [];
    
    % 更新分隔符位置
    idx0_new = find(offspring == 0);
    idx0_new = [0 idx0_new length(offspring)+1];
    
    % 插入到目标段的随机位置
    if toSeg <= length(idx0_new)-1
        toRange_new = (idx0_new(toSeg)+1):(idx0_new(toSeg+1)-1);
        if isempty(toRange_new)
            insertPos = idx0_new(toSeg) + 1;
        else
            insertPos = toRange_new(randi(length(toRange_new)));
        end
        offspring = [offspring(1:insertPos-1) pointVal offspring(insertPos:end)];
    else
        offspring = [offspring pointVal];
    end
    
    % 确保维度一致
    if length(offspring) ~= originalDim
        if length(offspring) < originalDim
            offspring = [offspring, zeros(1, originalDim - length(offspring))];
        else
            offspring = offspring(1:originalDim);
        end
    end
end
