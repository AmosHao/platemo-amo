function offspring = SimpleCrossover(parent1, parent2, m, varDim, n)
% 简化的交叉操作：保持0分隔符，确保非0值不重复（1-20每个值只出现一次）
%
% 输入:
%   parent1, parent2 - 父代解（包含0分隔符）
%   m - 无人机数量（段数）
%   varDim - 变量维度
%   n - 商家数量（用于确定必需的点：1到2*n）
%
% 输出:
%   offspring - 子代解（保证：m-1个0分隔符，1-2*n每个值只出现一次）

    % 确保输入维度正确
    if length(parent1) ~= varDim
        parent1 = parent1(1:min(varDim, length(parent1)));
        if length(parent1) < varDim
            parent1 = [parent1, zeros(1, varDim - length(parent1))];
        end
    end
    if length(parent2) ~= varDim
        parent2 = parent2(1:min(varDim, length(parent2)));
        if length(parent2) < varDim
            parent2 = [parent2, zeros(1, varDim - length(parent2))];
        end
    end
    
    % 找到0分隔符的位置（应该都是m-1个）
    idx0_1 = find(parent1 == 0);
    idx0_2 = find(parent2 == 0);
    
    % 如果0分隔符数量不对，使用parent1的0分隔符位置，并修正
    if length(idx0_1) ~= m - 1
        % 修正parent1的0分隔符
        if length(idx0_1) < m - 1
            nonZeroIdx = find(parent1 ~= 0);
            if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
                numNeeded = (m - 1) - length(idx0_1);
                insertPos = round(linspace(1, length(nonZeroIdx)-1, numNeeded+1));
                insertPos = insertPos(2:end);
                for i = length(insertPos):-1:1
                    pos = nonZeroIdx(insertPos(i));
                    parent1 = [parent1(1:pos), 0, parent1(pos+1:end)];
                end
            else
                parent1 = [parent1, zeros(1, (m-1) - length(idx0_1))];
            end
            idx0_1 = find(parent1 == 0);
        else
            idx0_1 = idx0_1(1:m-1);
        end
        % 调整维度
        if length(parent1) ~= varDim
            if length(parent1) < varDim
                nonZeroVals = parent1(parent1 ~= 0);
                if ~isempty(nonZeroVals)
                    fillValues = nonZeroVals(randi(length(nonZeroVals), 1, varDim - length(parent1)));
                    parent1 = [parent1, fillValues];
                end
            else
                parent1 = parent1(1:varDim);
            end
        end
        idx0_1 = find(parent1 == 0);
    end
    
    % 提取非0元素（去除重复，只保留唯一值）
    nonZero1 = parent1(parent1 ~= 0);
    nonZero2 = parent2(parent2 ~= 0);
    
    % 获取唯一值（去除重复）
    unique1 = unique(nonZero1, 'stable');
    unique2 = unique(nonZero2, 'stable');
    
    % 交叉策略：从parent1开始，用parent2的值替换一部分
    % 确保结果包含所有必需的点（1到2*n），且不重复
    allRequiredPoints = 1:(2*n);
    
    % 从parent1复制非0值（去除重复）
    offspring_nonZero = unique1;
    
    % 随机选择一部分位置，用parent2的值替换（但要确保不重复）
    if length(offspring_nonZero) > 0 && length(unique2) > 0
        numToReplace = randi([1, min(length(offspring_nonZero), length(unique2))]);
        replaceIdx = randperm(length(offspring_nonZero), numToReplace);
        
        % 从parent2选择值来替换（确保不重复）
        availableFrom2 = setdiff(unique2, offspring_nonZero);  % parent2中不在offspring的值
        if length(availableFrom2) > 0
            numToTake = min(numToReplace, length(availableFrom2));
            replaceVals = availableFrom2(randperm(length(availableFrom2), numToTake));
            
            % 替换
            for i = 1:min(length(replaceVals), length(replaceIdx))
                offspring_nonZero(replaceIdx(i)) = replaceVals(i);
            end
        end
    end
    
    % 检查是否包含所有必需的点
    missingPoints = setdiff(allRequiredPoints, offspring_nonZero);
    if ~isempty(missingPoints)
        % 补充缺失的点
        % 如果有重复值，先替换重复值
        [uniqueVals, ~, idx] = unique(offspring_nonZero);
        freq = histcounts(idx, 1:max(idx)+1);
        duplicateIdx = find(freq > 1);
        if ~isempty(duplicateIdx)
            % 替换重复值
            dupIdx = 1;
            for i = 1:length(offspring_nonZero)
                if ismember(offspring_nonZero(i), uniqueVals(duplicateIdx))
                    prevOccurrences = find(offspring_nonZero(1:i-1) == offspring_nonZero(i));
                    if ~isempty(prevOccurrences) && dupIdx <= length(missingPoints)
                        offspring_nonZero(i) = missingPoints(dupIdx);
                        dupIdx = dupIdx + 1;
                    end
                end
            end
        end
        
        % 如果还有缺失的点，添加到末尾（但要注意总长度）
        remainingMissing = setdiff(allRequiredPoints, offspring_nonZero);
        if ~isempty(remainingMissing)
            % 计算可以添加多少个点（考虑0分隔符）
            maxNonZero = varDim - (m - 1);  % 最多可以有多少个非0值
            numCanAdd = min(length(remainingMissing), maxNonZero - length(offspring_nonZero));
            if numCanAdd > 0
                offspring_nonZero = [offspring_nonZero, remainingMissing(1:numCanAdd)];
            end
        end
    end
    
    % 确保不重复（最终检查）
    [uniqueVals, ~, idx] = unique(offspring_nonZero);
    freq = histcounts(idx, 1:max(idx)+1);
    duplicateIdx = find(freq > 1);
    if ~isempty(duplicateIdx)
        % 如果有重复，只保留第一次出现的
        [~, firstOccurrence] = unique(offspring_nonZero, 'stable');
        offspring_nonZero = offspring_nonZero(sort(firstOccurrence));
    end
    
    % 限制长度（考虑0分隔符）
    maxNonZero = varDim - (m - 1);
    if length(offspring_nonZero) > maxNonZero
        offspring_nonZero = offspring_nonZero(1:maxNonZero);
    end
    
    % 重建解：保持0分隔符位置不变
    offspring = zeros(1, varDim);
    zeroPositions = idx0_1(1:min(m-1, length(idx0_1)));  % 确保不超过m-1个
    if length(zeroPositions) < m - 1
        % 如果0分隔符不足，在非0位置之间插入
        nonZeroPositions = setdiff(1:varDim, zeroPositions);
        if length(nonZeroPositions) > 1
            insertPos = round(linspace(1, length(nonZeroPositions)-1, (m-1)-length(zeroPositions)+1));
            insertPos = insertPos(2:end);
            for i = length(insertPos):-1:1
                pos = nonZeroPositions(insertPos(i));
                zeroPositions = [zeroPositions, pos];
            end
            zeroPositions = sort(zeroPositions);
        end
    end
    zeroPositions = zeroPositions(1:min(m-1, length(zeroPositions)));
    
    nonZeroPositions = setdiff(1:varDim, zeroPositions);
    
    % 将非0元素填入非0位置
    numNonZero = min(length(offspring_nonZero), length(nonZeroPositions));
    offspring(nonZeroPositions(1:numNonZero)) = offspring_nonZero(1:numNonZero);
    
    % 确保0分隔符位置正确
    offspring(zeroPositions) = 0;
    
    % 最终验证：确保0分隔符数量正确，非0值不重复
    zeroCount = sum(offspring == 0);
    if zeroCount ~= m - 1
        % 强制修正0分隔符
        if zeroCount < m - 1
            nonZeroIdx = find(offspring ~= 0);
            if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
                numNeeded = (m - 1) - zeroCount;
                insertPos = round(linspace(1, length(nonZeroIdx)-1, numNeeded+1));
                insertPos = insertPos(2:end);
                for i = length(insertPos):-1:1
                    pos = nonZeroIdx(insertPos(i));
                    offspring = [offspring(1:pos), 0, offspring(pos+1:end)];
                end
            end
        elseif zeroCount > m - 1
            zeroIdx = find(offspring == 0);
            offspring(zeroIdx(m:end)) = [];
        end
        % 调整维度
        if length(offspring) ~= varDim
            if length(offspring) < varDim
                nonZeroVals = offspring(offspring ~= 0);
                if ~isempty(nonZeroVals)
                    % 补充缺失的点，而不是重复已有的点
                    allPoints = 1:(2*n);
                    missing = setdiff(allPoints, nonZeroVals);
                    if ~isempty(missing)
                        numToAdd = min(length(missing), varDim - length(offspring));
                        offspring = [offspring, missing(1:numToAdd)];
                    else
                        % 如果没有缺失的点，用随机值填充（但确保不重复）
                        used = unique(nonZeroVals);
                        available = setdiff(1:(2*n), used);
                        if ~isempty(available)
                            numToAdd = min(length(available), varDim - length(offspring));
                            offspring = [offspring, available(randperm(length(available), numToAdd))];
                        end
                    end
                end
            else
                offspring = offspring(1:varDim);
            end
        end
    end
    
    % 最终验证：确保非0值不重复
    nonZeroVals = offspring(offspring ~= 0);
    if length(nonZeroVals) ~= length(unique(nonZeroVals))
        % 如果有重复，只保留第一次出现的
        [~, firstOccurrence] = unique(nonZeroVals, 'stable');
        nonZeroPositions = find(offspring ~= 0);
        offspring(nonZeroPositions) = 0;
        offspring(nonZeroPositions(sort(firstOccurrence))) = nonZeroVals(sort(firstOccurrence));
        
        % 补充缺失的点
        allPoints = 1:(2*n);
        missing = setdiff(allPoints, nonZeroVals(sort(firstOccurrence)));
        if ~isempty(missing)
            zeroPositions = find(offspring == 0);
            if length(zeroPositions) > m - 1
                % 用缺失的点替换多余的0
                numToReplace = min(length(missing), length(zeroPositions) - (m-1));
                replaceIdx = zeroPositions((m):(m-1+numToReplace));
                offspring(replaceIdx) = missing(1:numToReplace);
            end
        end
    end
end
