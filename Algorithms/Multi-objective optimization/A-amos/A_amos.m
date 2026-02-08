classdef A_amos < ALGORITHM
% <multi> <real/integer>
% A-amos: Adaptive Multi-objective Optimization for Order Assignment
% 基于PRMO改造的订单分配多目标优化算法
% 
% 参数说明:
% Kmax --- 12 --- 目标空间聚类最大数量（随种群增大而增加）
% BETA --- 0.8 --- 邻域选择平衡参数
% pc --- 0.9 --- 交叉概率（保持高交叉率以充分利用大种群）
% pm --- 0.2 --- 变异概率（适度增加以提高探索能力）
% maxgen --- 400 --- 最大迭代次数（增加迭代以充分利用大种群）

    methods
        function main(Algorithm, Problem)
           % 屏蔽复数绘图警告
           warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
           
           %% 参数设置
           % 调整参数以适应更大的种群：
           % - Kmax: 聚类数随种群增大而增加
           % - pc/pm: 保持较高的交叉概率，增加变异概率以提高多样性
           % - maxgen: 增加迭代次数以充分利用大种群
           [Kmax, BETA, pc, pm, maxgen] = Algorithm.ParameterSet(12, 0.8, 0.8, 0.3, 400);
           
           % 修复概率（降低以增加多样性，避免过度修复导致局部最优）
           repairProb = 0.6;  % 60%概率进行顺序约束修复
           completenessProb = 0.7;  % 70%概率进行完整性修复
           
           popSize = Problem.N; 
           objDim = Problem.M;
           varDim = Problem.D;
           n = Problem.n;  % 客户点数量（商家数量）
           m = Problem.m;  % 无人机数量
           
           %% 初始化种群
           Population = Problem.Initialization();
           pop = Population.decs;
           objVals = Population.objs;
           
           % 初始化聚类信息（在目标空间）
           clustTag = (1:popSize)';
           clustName = (1:popSize)';
           centroid = objVals;  % 聚类中心在目标空间
           
           %% 主优化循环
           gen = 0;
           while Algorithm.NotTerminated(Population)
               gen = gen + 1;
               
               %% 目标空间聚类管理
               auxPop = [];
               globalClust = [];  % 获取全局子种群
               nClust = size(clustName, 1);  % 聚类数量
               neigSet = cell(nClust, 1);
               
               % 确保pop和objVals不为空
               if isempty(pop) || isempty(objVals)
                   continue;  % 如果种群为空，跳过本次迭代
               end
               
               % 为每个聚类构建邻域集合（基于目标空间距离）
               for i = 1:nClust
                   if i > length(clustName)
                       continue;  % 跳过无效索引
                   end
                   mark = clustTag == clustName(i);  % 找到当前属于第i类的个体编号
                   if sum(mark) > 0
                       neigPop = pop(mark, :);
                       neigObj = objVals(mark, :);
                       
                       % 从每个聚类随机选择一个个体加入到全局子种群
                       pos = randsample(sum(mark), 1);
                       globalClust = [globalClust; neigPop(pos, :)];
                       
                       if sum(mark) > 2
                           neigSet{i, 1} = {neigPop, neigObj};  % 存储邻域解和目标值
                       end
                   end
               end
               globalSize = size(globalClust, 1);
               
               %% 个体进化
               for i = 1:popSize
                   if i > size(pop, 1) || i > size(objVals, 1) || i > length(clustTag)
                       continue;  % 跳过无效索引
                   end
                   currentSol = pop(i, :);
                   currentObj = objVals(i, :);
                   currentTag = clustTag(i);
                   
                   % 确保currentSol维度正确
                   if length(currentSol) ~= varDim
                       if length(currentSol) < varDim
                           currentSol = [currentSol, zeros(1, varDim - length(currentSol))];
                       else
                           currentSol = currentSol(1:varDim);
                       end
                   end
                   
                   % 获取当前个体的邻域
                   neighborhood = [];
                   neighborhoodObj = [];
                   clustIdx = find(clustName == currentTag, 1);
                   if ~isempty(clustIdx) && clustIdx <= length(neigSet) && ~isempty(neigSet{clustIdx, 1})
                       temp = neigSet{clustIdx, 1};
                       neighborhood = temp{1};
                       neighborhoodObj = temp{2};
                       % 从邻域中删除当前个体
                       if ~isempty(neighborhood)
                           [~, loc] = ismember(currentSol, neighborhood, 'rows');
                           if loc > 0 && loc <= size(neighborhood, 1)
                               neighborhood(loc, :) = [];
                               neighborhoodObj(loc, :) = [];
                           end
                       end
                   end
                   neigSize = size(neighborhood, 1);
                   
                   %% 父代选择策略：以概率选择从全局还是当前解的聚类中选择一个个体
                   % 聚类的作用：选择同聚类的个体进行局部搜索，或选择全局个体进行探索
                   localSearchProb = BETA * gen / maxgen;  % 局部搜索概率（随代数增加）
                   
                   % 以概率选择从当前解的聚类还是全局选择
                   if rand < localSearchProb && neigSize >= 1
                       % 局部搜索：从当前解的聚类（邻域）中选择一个个体
                       idx = randsample(neigSize, 1);
                       parent = neighborhood(idx, :);
                   elseif globalSize >= 1
                       % 全局搜索：从全局聚类中选择一个个体
                       idx = randsample(globalSize, 1);
                       parent = globalClust(idx, :);
                   else
                       % 备用：使用当前解
                       parent = currentSol;
                   end
                   
                   %% 对选中的父代应用变异操作生成子代
                   % 始终应用变异操作，确保生成新的子代
                   trialSol = OrderPreservingMutation(parent, m, varDim, n);
                   
                   % 确保是行向量
                   if size(trialSol, 1) ~= 1
                       trialSol = reshape(trialSol, 1, []);
                   end
                   
                   % 添加到auxPop
                   if isempty(auxPop)
                       auxPop = trialSol;
                   else
                       auxPop = [auxPop; trialSol];
                   end
               end
               
               %% 环境选择
               selectedSize = Problem.N;
               [pop, objVals, clustTag, clustName, centroid, Population] = ...
                   EnvironmentalSelection(auxPop, pop, objVals, clustTag, clustName, ...
                   centroid, selectedSize, objDim, varDim, Kmax, Problem);
               
               popSize = size(pop, 1);
               
               % 更新Population对象
               if popSize > 0 && length(Population) >= popSize
                   for i = 1:popSize
                       if i <= length(Population)
                           Population(i).add_clustTag = clustTag(i, :);
                           Population(i).add_clustName = clustName(:, :);
                       end
                   end
               end
           end
        end
    end
end

function completeSol = EnsureSolutionCompleteness(sol, n, m, varDim)
% 确保解包含所有必需的点（1到2*n）
% **关键修复**：确保0分隔符数量正确（m-1个），并在调整维度时优先保留0分隔符

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
