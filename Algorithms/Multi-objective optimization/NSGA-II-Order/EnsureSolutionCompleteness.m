function completeSol = EnsureSolutionCompleteness(sol, n, m, varDim)
% 确保解包含所有必需的点（1到2*n）
% 这是解决无效解过多的关键：不完整的解会导致目标值异常

    completeSol = sol;
    
    % 必需的点：1到2*n（商家1-10，客户点11-20）
    allRequiredPoints = 1:(2*n);
    
    % 找到所有非0的点
    nonZeroVals = completeSol(completeSol ~= 0);
    uniqueVals = unique(nonZeroVals);
    
    % 检查是否缺少必需的点
    missingVals = setdiff(allRequiredPoints, uniqueVals);
    
    if ~isempty(missingVals)
        % 如果缺少点，需要补充
        % 策略：替换重复的点，或插入到0分隔符附近
        
        % 找出重复的点
        if ~isempty(nonZeroVals)
            [counts, ~, idx] = unique(nonZeroVals);
            freq = histcounts(idx, 1:max(idx)+1);
            duplicateVals = nonZeroVals(freq > 1);
        else
            duplicateVals = [];
        end
        
        % 先替换重复的点
        dupIdx = 1;
        for i = 1:length(completeSol)
            if completeSol(i) ~= 0 && ismember(completeSol(i), duplicateVals)
                % 检查是否已经出现过
                prevOccurrences = find(completeSol(1:i-1) == completeSol(i));
                if ~isempty(prevOccurrences) && dupIdx <= length(missingVals)
                    completeSol(i) = missingVals(dupIdx);
                    dupIdx = dupIdx + 1;
                end
            end
        end
        
        % 如果还有缺失的点，插入到0分隔符附近
        remainingMissing = setdiff(allRequiredPoints, unique(completeSol(completeSol ~= 0)));
        if ~isempty(remainingMissing)
            zeroIdx = find(completeSol == 0);
            if ~isempty(zeroIdx) && zeroIdx(1) > 1
                % 在第一个0之前插入
                completeSol = [completeSol(1:zeroIdx(1)-1), remainingMissing, completeSol(zeroIdx(1):end)];
            else
                % 在末尾插入
                completeSol = [completeSol, remainingMissing];
            end
        end
    end
    
    % **关键修复**：先确保0分隔符数量正确（m-1个），再调整维度
    % 这样可以避免在调整维度时丢失0分隔符
    zeroCount = sum(completeSol == 0);
    if zeroCount < m - 1
        % 如果0分隔符不足，需要插入
        % 策略：在非0元素之间均匀插入0分隔符
        nonZeroIdx = find(completeSol ~= 0);
        if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
            % 计算需要插入的0的数量
            numZerosNeeded = (m - 1) - zeroCount;
            % 在非0元素之间均匀插入0
            insertPositions = round(linspace(1, length(nonZeroIdx)-1, numZerosNeeded+1));
            insertPositions = insertPositions(2:end);  % 去掉第一个位置
            % 从后往前插入，避免索引变化
            for i = length(insertPositions):-1:1
                pos = nonZeroIdx(insertPositions(i));
                completeSol = [completeSol(1:pos), 0, completeSol(pos+1:end)];
            end
        else
            % 如果没有非0元素或只有一个，直接在末尾添加0
            completeSol = [completeSol, zeros(1, (m-1) - zeroCount)];
        end
    elseif zeroCount > m - 1
        % 如果0分隔符过多，只保留前m-1个
        zeroIdx = find(completeSol == 0);
        if length(zeroIdx) > m - 1
            % 删除多余的0，但要保留前m-1个
            toRemove = zeroIdx(m:end);
            completeSol(toRemove) = [];
        end
    end
    
    % 确保维度正确（在确保0分隔符之后）
    % **关键**：调整维度时，要确保0分隔符不被截断
    currentZeroCount = sum(completeSol == 0);
    if length(completeSol) ~= varDim
        if length(completeSol) < varDim
            % 如果长度不足，补充0（但要确保0分隔符数量不超过m-1）
            remainingZeros = (m-1) - currentZeroCount;
            if remainingZeros > 0
                % 先补充0分隔符
                completeSol = [completeSol, zeros(1, remainingZeros)];
            end
            % 再补充到目标维度（用非0值填充，但这里我们用0，因为varDim应该已经包含了0分隔符）
            % 实际上，varDim应该等于 2*n + (m-1)，所以如果长度不足，应该是缺少非0元素
            % 但为了安全，我们只补充到varDim，不再添加额外的0
            if length(completeSol) < varDim
                completeSol = [completeSol, zeros(1, varDim - length(completeSol))];
            end
        else
            % 如果长度超过varDim，需要截断
            % **关键**：截断时要确保保留m-1个0分隔符
            zeroIdx = find(completeSol == 0);
            if length(zeroIdx) >= m - 1
                % 保留前m-1个0分隔符，截断其他部分
                lastZeroIdx = zeroIdx(m-1);
                % 计算需要保留的非0元素数量
                nonZeroNeeded = varDim - (m - 1);
                nonZeroCount = sum(completeSol(1:lastZeroIdx) ~= 0);
                if nonZeroCount <= nonZeroNeeded
                    % 保留到最后一个0分隔符之后，再补充一些非0元素
                    completeSol = completeSol(1:min(varDim, length(completeSol)));
                else
                    % 需要截断一些非0元素
                    completeSol = completeSol(1:varDim);
                end
            else
                % 如果0分隔符不足，先补充0，再截断
                completeSol = completeSol(1:varDim);
            end
        end
    end
    
    % 最终验证：确保0分隔符数量正确
    finalZeroCount = sum(completeSol == 0);
    if finalZeroCount ~= m - 1
        % 如果还是不对，强制修正
        zeroIdx = find(completeSol == 0);
        if finalZeroCount < m - 1
            % 补充0分隔符
            numNeeded = (m - 1) - finalZeroCount;
            % 在非0元素之间插入
            nonZeroIdx = find(completeSol ~= 0);
            if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
                insertPos = round(linspace(1, length(nonZeroIdx)-1, numNeeded+1));
                insertPos = insertPos(2:end);
                for i = length(insertPos):-1:1
                    pos = nonZeroIdx(insertPos(i));
                    completeSol = [completeSol(1:pos), 0, completeSol(pos+1:end)];
                end
            else
                completeSol = [completeSol, zeros(1, numNeeded)];
            end
        elseif finalZeroCount > m - 1
            % 删除多余的0
            zeroIdx = find(completeSol == 0);
            completeSol(zeroIdx(m:end)) = [];
        end
        % 再次调整维度
        if length(completeSol) ~= varDim
            if length(completeSol) < varDim
                completeSol = [completeSol, zeros(1, varDim - length(completeSol))];
            else
                completeSol = completeSol(1:varDim);
            end
        end
    end
end
