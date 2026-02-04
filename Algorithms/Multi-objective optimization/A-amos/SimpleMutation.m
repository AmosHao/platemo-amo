function offspring = SimpleMutation(individual, m, varDim, n)
% 简化的变异操作：只操作非0元素，保持0分隔符位置和数量不变，确保非0值不重复
%
% 输入:
%   individual - 待变异的个体
%   m - 无人机数量（段数）
%   varDim - 变量维度
%   n - 商家数量（用于确定必需的点：1到2*n）
%
% 输出:
%   offspring - 变异后的个体（保证：m-1个0分隔符，1-2*n每个值只出现一次）

    % 确保输入维度正确
    if length(individual) ~= varDim
        individual = individual(1:min(varDim, length(individual)));
        if length(individual) < varDim
            individual = [individual, zeros(1, varDim - length(individual))];
        end
    end
    
    % 找到0分隔符的位置
    idx0 = find(individual == 0);
    
    % 如果0分隔符数量不对，先修正
    if length(idx0) ~= m - 1
        if length(idx0) < m - 1
            nonZeroIdx = find(individual ~= 0);
            if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
                numNeeded = (m - 1) - length(idx0);
                insertPos = round(linspace(1, length(nonZeroIdx)-1, numNeeded+1));
                insertPos = insertPos(2:end);
                for i = length(insertPos):-1:1
                    pos = nonZeroIdx(insertPos(i));
                    individual = [individual(1:pos), 0, individual(pos+1:end)];
                end
            else
                individual = [individual, zeros(1, (m-1) - length(idx0))];
            end
            idx0 = find(individual == 0);
        else
            idx0 = idx0(1:m-1);
        end
        % 调整维度
        if length(individual) ~= varDim
            if length(individual) < varDim
                nonZeroVals = individual(individual ~= 0);
                allPoints = 1:(2*n);
                missing = setdiff(allPoints, nonZeroVals);
                if ~isempty(missing)
                    numToAdd = min(length(missing), varDim - length(individual));
                    individual = [individual, missing(1:numToAdd)];
                end
            else
                individual = individual(1:varDim);
            end
        end
        idx0 = find(individual == 0);
    end
    
    % 提取非0元素
    nonZeroIdx = find(individual ~= 0);
    nonZeroVals = individual(nonZeroIdx);
    
    if length(nonZeroVals) < 2
        offspring = individual;
        return;
    end
    
    % 简单的交换变异：随机选择两个非0位置交换
    swapIdx = randperm(length(nonZeroIdx), 2);
    
    offspring = individual;
    temp = offspring(nonZeroIdx(swapIdx(1)));
    offspring(nonZeroIdx(swapIdx(1))) = offspring(nonZeroIdx(swapIdx(2)));
    offspring(nonZeroIdx(swapIdx(2))) = temp;
    
    % 验证：0分隔符应该保持不变
    if sum(offspring == 0) ~= m - 1
        % 如果0分隔符数量改变，恢复
        offspring = individual;
    end
    
    % 验证：确保非0值不重复
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
