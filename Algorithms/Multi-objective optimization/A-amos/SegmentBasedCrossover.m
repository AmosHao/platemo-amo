function offspring = SegmentBasedCrossover(parent1, parent2, m)
% 基于段的交叉操作
% 保留0分隔符结构，在段级别进行交叉
%
% 输入:
%   parent1, parent2 - 父代解（包含0分隔符）
%   m - 无人机数量（段数）
%
% 输出:
%   offspring - 子代解

    % 保存原始维度
    originalDim1 = length(parent1);
    originalDim2 = length(parent2);
    expectedDim = max(originalDim1, originalDim2);  % 使用较大的维度
    
    % 找到0分隔符的位置
    idx0_1 = find(parent1 == 0);
    idx0_2 = find(parent2 == 0);
    
    % 确保有 m-1 个0分隔符
    if length(idx0_1) < m - 1
        % 如果0的数量不足，在末尾补充0
        numZerosNeeded = (m - 1) - length(idx0_1);
        parent1 = [parent1, zeros(1, numZerosNeeded)];
        idx0_1 = find(parent1 == 0);
    elseif length(idx0_1) > m - 1
        % 如果0的数量过多，只保留前m-1个
        idx0_1 = idx0_1(1:m-1);
    end
    
    if length(idx0_2) < m - 1
        % 如果0的数量不足，在末尾补充0
        numZerosNeeded = (m - 1) - length(idx0_2);
        parent2 = [parent2, zeros(1, numZerosNeeded)];
        idx0_2 = find(parent2 == 0);
    elseif length(idx0_2) > m - 1
        % 如果0的数量过多，只保留前m-1个
        idx0_2 = idx0_2(1:m-1);
    end
    
    % 将解分段
    segments1 = cell(m, 1);
    segments2 = cell(m, 1);
    
    idx0_1 = [0 idx0_1 length(parent1)+1];
    idx0_2 = [0 idx0_2 length(parent2)+1];
    
    % 确保idx0的长度足够
    for k = 1:m
        if k < length(idx0_1)
            segments1{k} = parent1(idx0_1(k)+1 : idx0_1(k+1)-1);
        else
            segments1{k} = parent1(idx0_1(k)+1 : idx0_1(end)-1);
        end
        
        if k < length(idx0_2)
            segments2{k} = parent2(idx0_2(k)+1 : idx0_2(k+1)-1);
        else
            segments2{k} = parent2(idx0_2(k)+1 : idx0_2(end)-1);
        end
    end
    
    % 随机选择要交换的段
    if m > 1
        numSegsToSwap = randi([1, min(m-1, m)]);  % 至少交换1段，最多交换m-1段
        numSegsToSwap = min(numSegsToSwap, m);  % 确保不超过m
        segsToSwap = randperm(m, numSegsToSwap);
    else
        % 如果只有1段，不交换
        segsToSwap = [];
    end
    
    % 创建子代
    offspring = [];
    for k = 1:m
        if ismember(k, segsToSwap)
            % 从parent2取段
            offspring = [offspring segments2{k}];
        else
            % 从parent1取段
            offspring = [offspring segments1{k}];
        end
        
        % 添加0分隔符（除了最后一段）
        if k < m
            offspring = [offspring 0];
        end
    end
    
    % 修复重复和缺失的点
    offspring = FixSequence(offspring, m);
    
    % 确保维度与父代一致（使用原始维度）
    if length(offspring) ~= expectedDim
        if length(offspring) < expectedDim
            offspring = [offspring, zeros(1, expectedDim - length(offspring))];
        else
            offspring = offspring(1:expectedDim);
        end
    end
end

function fixed = FixSequence(seq, m)
% 修复序列中的重复和缺失点
% 保持0分隔符位置不变
% **关键修复**：强制确保包含所有必需的点（1到2*n，其中n是商家数量）
% 这是解决局部最优的关键：不完整的解会导致目标值异常

    % 找到所有非0的点
    nonZeroIdx = find(seq ~= 0);
    nonZeroVals = seq(nonZeroIdx);
    
    if isempty(nonZeroVals)
        fixed = seq;
        return;
    end
    
    % **关键修复**：强制要求包含所有1到2*n的点（n=10，所以是1-20）
    % 假设n=10（商家数量），所以必需的点是1-20
    % 如果不知道n，使用max(nonZeroVals)作为上限，但通常应该是2*n
    allRequiredPoints = 1:max(20, max(nonZeroVals));  % 至少1-20，或更大范围
    
    % 找出重复和缺失的点
    uniqueVals = unique(nonZeroVals);
    missingVals = setdiff(allRequiredPoints, uniqueVals);
    
    % 找出重复的点
    [counts, ~, idx] = unique(nonZeroVals);
    freq = histcounts(idx, 1:max(idx)+1);
    duplicateVals = nonZeroVals(freq > 1);
    
    % 第一步：替换重复的点
    fixed = seq;
    dupIdx = 1;
    for i = 1:length(nonZeroIdx)
        if ismember(nonZeroVals(i), duplicateVals)
            % 检查这个值是否已经出现过
            prevOccurrences = find(fixed(1:nonZeroIdx(i)-1) == nonZeroVals(i));
            if ~isempty(prevOccurrences)
                % 这是一个重复值，需要替换为缺失的点
                if dupIdx <= length(missingVals)
                    fixed(nonZeroIdx(i)) = missingVals(dupIdx);
                    dupIdx = dupIdx + 1;
                else
                    % 如果没有缺失值，随机替换为未使用的值
                    usedVals = unique(fixed(fixed ~= 0));
                    unusedVals = setdiff(allRequiredPoints, usedVals);
                    if ~isempty(unusedVals)
                        fixed(nonZeroIdx(i)) = unusedVals(randi(length(unusedVals)));
                    end
                end
            end
        end
    end
    
    % 第二步：补充缺失的点（关键修复）
    % 更新缺失的点列表（因为可能已经替换了一些重复点）
    currentNonZeroVals = fixed(fixed ~= 0);
    currentUniqueVals = unique(currentNonZeroVals);
    remainingMissingVals = setdiff(allRequiredPoints, currentUniqueVals);
    
    if ~isempty(remainingMissingVals)
        % 找到0分隔符的位置
        zeroIdx = find(fixed == 0);
        
        if ~isempty(zeroIdx)
            % 在0分隔符附近插入缺失的点
            % 策略：尽量均匀分布到各个段
            numSegs = length(zeroIdx) + 1;
            pointsPerSeg = ceil(length(remainingMissingVals) / numSegs);
            
            insertIdx = 1;
            for seg = 1:numSegs
                if insertIdx > length(remainingMissingVals)
                    break;
                end
                
                % 确定插入位置
                if seg == 1
                    % 第一段：在第一个0之前
                    if zeroIdx(1) > 1
                        insertPos = zeroIdx(1);
                        numToInsert = min(pointsPerSeg, length(remainingMissingVals) - insertIdx + 1);
                        fixed = [fixed(1:insertPos-1), remainingMissingVals(insertIdx:insertIdx+numToInsert-1), fixed(insertPos:end)];
                        insertIdx = insertIdx + numToInsert;
                        zeroIdx = find(fixed == 0);  % 更新0的位置
                    end
                elseif seg <= length(zeroIdx)
                    % 中间段：在0之后
                    insertPos = zeroIdx(seg);
                    numToInsert = min(pointsPerSeg, length(remainingMissingVals) - insertIdx + 1);
                    fixed = [fixed(1:insertPos), remainingMissingVals(insertIdx:insertIdx+numToInsert-1), fixed(insertPos+1:end)];
                    insertIdx = insertIdx + numToInsert;
                    zeroIdx = find(fixed == 0);  % 更新0的位置
                else
                    % 最后一段：在末尾
                    numToInsert = length(remainingMissingVals) - insertIdx + 1;
                    fixed = [fixed, remainingMissingVals(insertIdx:end)];
                    break;
                end
            end
        else
            % 如果没有0分隔符，在末尾插入
            fixed = [fixed, remainingMissingVals];
        end
    end
    
    % 确保0分隔符数量正确（m-1个）
    zeroCount = sum(fixed == 0);
    if zeroCount < m - 1
        % 补充0
        fixed = [fixed, zeros(1, (m-1) - zeroCount)];
    elseif zeroCount > m - 1
        % 移除多余的0（保留前m-1个）
        zeroIdx = find(fixed == 0);
        if length(zeroIdx) > m - 1
            fixed(zeroIdx(m:end)) = [];
        end
    end
end
