function offspring = SegmentBasedCrossover(parent1, parent2, n, m, varDim)
% 基于三种策略的交叉入口：随机选用位置追踪 / 对优先级 / 顺序保持 之一
% 保证：0 个数为 m-1，每段内商家-客户对完整，全局 1..2n 各出现一次，段内商家在客户前
%
% 输入:
%   parent1, parent2 - 父代解（包含0分隔符，长度宜为 varDim）
%   n - 商家数量（客户点数也为 n，编码 1..n 商家、n+1..2n 客户）
%   m - 无人机数量（段数，0 的个数为 m-1）
%   varDim - 决策变量维度（应为 2*n + (m-1)）
%
% 输出:
%   offspring - 子代解（长度 varDim，恰好 m-1 个 0，完整 1..2n，段内不分割对）

    % 统一父代长度为 varDim（便于后续交叉与修复）
    if length(parent1) ~= varDim
        if length(parent1) < varDim
            parent1 = [parent1, zeros(1, varDim - length(parent1))];
        else
            parent1 = parent1(1:varDim);
        end
    end
    if length(parent2) ~= varDim
        if length(parent2) < varDim
            parent2 = [parent2, zeros(1, varDim - length(parent2))];
        else
            parent2 = parent2(1:varDim);
        end
    end

    % 自适应选择交叉方式：停滞时偏向探索性更强的策略
    % 使用轮盘赌选择，但根据情况调整权重
    persistent crossoverWeights crossoverCounts;
    if isempty(crossoverWeights)
        crossoverWeights = [1, 1, 1];  % 初始权重相等
        crossoverCounts = [0, 0, 0];
    end
    
    % 如果传入停滞信息（通过varDim的符号位或额外参数），调整权重
    % 这里简化：随机选择但增加探索性
    if rand < 0.4
        % 40%概率使用位置追踪（探索性强）
        offspring = PositionTrackingCrossover(parent1, parent2, n, m, varDim);
        crossoverCounts(1) = crossoverCounts(1) + 1;
    elseif rand < 0.7
        % 30%概率使用对优先级（平衡）
        offspring = PairPriorityCrossover(parent1, parent2, n, m, varDim);
        crossoverCounts(2) = crossoverCounts(2) + 1;
    else
        % 30%概率使用顺序保持（保守）
        offspring = OrderPreservingCrossover(parent1, parent2, n, m, varDim);
        crossoverCounts(3) = crossoverCounts(3) + 1;
    end

    % 最终保证长度与 0 个数
    if length(offspring) ~= varDim
        if length(offspring) < varDim
            offspring = [offspring, zeros(1, varDim - length(offspring))];
        else
            offspring = offspring(1:varDim);
        end
    end
    % 确保首尾不为 0（0 仅作段间分隔，首尾必须是商家或客户点）
    offspring = SegmentBasedCrossover_ensureNoZeroAtEnds(offspring);
end

function offspring = PositionTrackingCrossover(parent1, parent2, n, m, varDim)
% 位置追踪交叉：追踪商家-客户对的位置关系，确保交叉后保持顺序约束
%
% 输入:
%   parent1, parent2 - 父代个体（行向量）
%   n - 商家数量（客户点数量）
%   m - 无人机数量
%   varDim - 变量维度
%
% 输出:
%   offspring - 交叉后的子代个体

    % 1. 验证父代有效性
    if ~SegmentBasedCrossover_isValid(parent1, n) || ~SegmentBasedCrossover_isValid(parent2, n)
        offspring = parent1;
        return;
    end
    
    len = varDim;  % 决策变量固定长度
    if length(parent1) < len
        parent1 = [parent1, zeros(1, len - length(parent1))];
    else
        parent1 = parent1(1:len);
    end
    if length(parent2) < len
        parent2 = [parent2, zeros(1, len - length(parent2))];
    else
        parent2 = parent2(1:len);
    end
    
    % 2. 构建位置映射
    [pos_map1, valid1] = SegmentBasedCrossover_buildPositionMap(parent1, n);
    [pos_map2, valid2] = SegmentBasedCrossover_buildPositionMap(parent2, n);
    
    if ~valid1 || ~valid2
        offspring = parent1;
        return;
    end
    
    % 3. 选择安全的交叉点
    valid_crossover_points = [];
    for k = 1:len-1
        if SegmentBasedCrossover_isSafeCrossoverPoint(pos_map1, k, n) && SegmentBasedCrossover_isSafeCrossoverPoint(pos_map2, k, n)
            valid_crossover_points = [valid_crossover_points, k];
        end
    end
    
    if isempty(valid_crossover_points)
        offspring = parent1;
        return;
    end
    
    k = valid_crossover_points(randi(length(valid_crossover_points)));
    
    % 4. 生成子代
    offspring = zeros(1, len);
    
    % 前k位：复制父代1的前k位
    offspring(1:k) = parent1(1:k);
    
    % 后len-k位：从父代2中填充未使用的基因
    used = unique(offspring(1:k));
    fill_idx = k+1;
    
    for i = 1:len
        gene = parent2(i);
        if ~ismember(gene, used)
            if gene == 0
                % 0分隔符：直接添加
                offspring(fill_idx) = gene;
                used = [used, gene];
                fill_idx = fill_idx + 1;
            elseif gene > n && gene <= 2*n
                % 客户点：检查对应商家是否已添加
                merchant = gene - n;
                if ismember(merchant, used)
                    offspring(fill_idx) = gene;
                    used = [used, gene];
                    fill_idx = fill_idx + 1;
                end
            else
                % 商家：直接添加
                offspring(fill_idx) = gene;
                used = [used, gene];
                fill_idx = fill_idx + 1;
            end
        end
        if fill_idx > len
            break;
        end
    end
    
    % 5. 修复：确保包含所有必需的点和正确数量的0分隔符
    offspring = SegmentBasedCrossover_repairSolution(offspring, n, m, varDim);
    
    % 6. 验证并修复顺序约束
    if ~SegmentBasedCrossover_isValid(offspring, n)
        offspring = SegmentBasedCrossover_repairOrder(offspring, n);
    end
    % 7. 修复段内完整对：每段内不分割商家-客户对
    offspring = SegmentBasedCrossover_repairSegments(offspring, n, m);
end

function offspring = PairPriorityCrossover(parent1, parent2, n, m, varDim)
% 对优先级交叉：基于商家优先级构建子代，确保顺序约束
%
% 输入:
%   parent1, parent2 - 父代个体（行向量）
%   n - 商家数量（客户点数量）
%   m - 无人机数量
%   varDim - 变量维度
%
% 输出:
%   offspring - 交叉后的子代个体

    % 1. 验证父代有效性
    if ~SegmentBasedCrossover_isValid(parent1, n) || ~SegmentBasedCrossover_isValid(parent2, n)
        offspring = parent1;
        return;
    end
    
    len = length(parent1);
    
    % 2. 计算商家优先级（位置越前，优先级越高）
    priority1 = SegmentBasedCrossover_calculatePriority(parent1, n);
    priority2 = SegmentBasedCrossover_calculatePriority(parent2, n);
    
    % 3. 合并优先级（取最小值，优先级更高）
    merged_priority = min(priority1, priority2);
    
    % 4. 按优先级排序商家
    [~, sorted_merchants] = sort(merged_priority);
    
    % 5. 构建子代（长度 varDim，0 的个数由 repairSolution 统一为 m-1）
    offspring = zeros(1, varDim);
    len = varDim;
    current_pos = 1;
    
    % 5.1 放置商家（按优先级顺序）
    for i = 1:length(sorted_merchants)
        merchant = sorted_merchants(i);
        if current_pos <= len
            offspring(current_pos) = merchant;
            current_pos = current_pos + 1;
        end
    end
    
    % 5.2 放置客户点（在对应商家之后的第一个 0 位）
    for i = 1:length(sorted_merchants)
        merchant = sorted_merchants(i);
        customer = merchant + n;
        merchant_pos = find(offspring == merchant, 1);
        if ~isempty(merchant_pos) && merchant_pos < len
            for j = merchant_pos+1:len
                if offspring(j) == 0
                    offspring(j) = customer;
                    break;
                end
            end
        end
    end
    
    % 6. 修复：确保包含所有必需的点、m-1 个 0、长度 varDim
    offspring = SegmentBasedCrossover_repairSolution(offspring, n, m, varDim);
    
    % 7. 验证并修复顺序约束
    if ~SegmentBasedCrossover_isValid(offspring, n)
        offspring = SegmentBasedCrossover_repairOrder(offspring, n);
    end
    % 8. 修复段内完整对
    offspring = SegmentBasedCrossover_repairSegments(offspring, n, m);
end

function offspring = OrderPreservingCrossover(parent1, parent2, n, m, varDim)
% 顺序保持交叉：基于商家顺序构建子代，确保顺序约束
%
% 输入:
%   parent1, parent2 - 父代个体（行向量）
%   n - 商家数量（客户点数量）
%   m - 无人机数量
%   varDim - 变量维度
%
% 输出:
%   offspring - 交叉后的子代个体

    % 1. 验证父代有效性
    if ~SegmentBasedCrossover_isValid(parent1, n) || ~SegmentBasedCrossover_isValid(parent2, n)
        offspring = parent1;
        return;
    end
    
    len = length(parent1);
    
    % 2. 提取商家顺序
    order1 = SegmentBasedCrossover_extractMerchantOrder(parent1, n);
    order2 = SegmentBasedCrossover_extractMerchantOrder(parent2, n);
    
    if isempty(order1) || isempty(order2)
        offspring = parent1;
        return;
    end
    
    % 3. 融合商家顺序：融合点前取自 order1，融合点后取自 order2，再补全两边未出现的商家
    fusion_point = randi(length(order1));
    fused_order = [order1(1:fusion_point), order2(fusion_point+1:end)];
    fused_order = unique(fused_order, 'stable');
    % 补全：确保 1..n 的商家都出现（从两个父代顺序中补缺失）
    all_merchants = unique([order1, order2]);
    fused_order = [fused_order, setdiff(all_merchants, fused_order)];
    
    % 4. 构建子代（使用 varDim 保证长度，0 的个数由 repairSolution 统一）
    len = varDim;
    offspring = zeros(1, len);
    current_pos = 1;
    
    % 4.1 放置商家（按融合顺序）
    for i = 1:length(fused_order)
        merchant = fused_order(i);
        if current_pos <= len
            offspring(current_pos) = merchant;
            current_pos = current_pos + 1;
        end
    end
    
    % 4.2 放置客户点（在对应商家之后的第一个 0 位）
    for i = 1:length(fused_order)
        merchant = fused_order(i);
        customer = merchant + n;
        merchant_pos = find(offspring == merchant, 1);
        if ~isempty(merchant_pos) && merchant_pos < len
            for j = merchant_pos+1:len
                if offspring(j) == 0
                    offspring(j) = customer;
                    break;
                end
            end
        end
    end
    
    % 5. 修复：确保包含所有必需的点、m-1 个 0、长度 varDim
    offspring = SegmentBasedCrossover_repairSolution(offspring, n, m, varDim);
    
    % 6. 验证并修复顺序约束
    if ~SegmentBasedCrossover_isValid(offspring, n)
        offspring = SegmentBasedCrossover_repairOrder(offspring, n);
    end
    % 7. 修复段内完整对
    offspring = SegmentBasedCrossover_repairSegments(offspring, n, m);
end

function valid = SegmentBasedCrossover_isValid(individual, n)
% 辅助函数：验证解的有效性（商家在对应客户点之前）
    valid = true;
    for s = 1:n
        c = s + n;
        pos_s = find(individual == s, 1);
        pos_c = find(individual == c, 1);
        if isempty(pos_s) || isempty(pos_c) || pos_s > pos_c
            valid = false;
            return;
        end
    end
end

function [pos_map, valid] = SegmentBasedCrossover_buildPositionMap(individual, n)
% 辅助函数：构建位置映射
    pos_map = zeros(2*n, 1);
    valid = true;
    
    for i = 1:length(individual)
        gene = individual(i);
        if gene >= 1 && gene <= 2*n
            pos_map(gene) = i;
        end
    end
    
    % 验证所有商家在对应客户点之前
    for s = 1:n
        c = s + n;
        if pos_map(s) == 0 || pos_map(c) == 0 || pos_map(s) > pos_map(c)
            valid = false;
            return;
        end
    end
end

function repaired = SegmentBasedCrossover_repairSegmentSizes(seq, m)
% 辅助函数：确保每个 0 分隔的段至少 2 个基因（至少一对商家-客户），避免出现 "0 18" 这种末段只有一点
    pos0 = find(seq == 0);
    if isempty(pos0) || length(pos0) ~= m - 1
        repaired = seq;
        return;
    end
    L = length(seq);
    % 各段长度：段1 = pos0(1)-1, 段k = pos0(k)-pos0(k-1)-1, 段m = L - pos0(m-1)
    segLen = [pos0(1)-1, pos0(2:end)-pos0(1:end-1)-1, L - pos0(end)];
    maxIter = L;  % 防止死循环
    iter = 0;
    while any(segLen < 2) && iter < maxIter
        iter = iter + 1;
        if segLen(1) < 2
            % 第一段太短：把第一个 0 向右移（与右边基因交换）
            p = pos0(1);
            if p < L
                seq([p, p+1]) = seq([p+1, p]);
            end
        elseif segLen(m) < 2
            % 最后一段太短：把最后一个 0 向左移（与左边基因交换）
            p = pos0(end);
            if p > 1
                seq([p-1, p]) = seq([p, p-1]);
            end
        else
            % 中间某段太短
            k = find(segLen < 2, 1);
            pLeft = pos0(k-1);   % 该段左边界 0
            pRight = pos0(k);    % 该段右边界 0
            % 若左边段更长，把左 0 右移；否则把右 0 左移
            if k > 1 && segLen(k-1) >= 2
                if pLeft < L
                    seq([pLeft, pLeft+1]) = seq([pLeft+1, pLeft]);
                end
            elseif k < m && segLen(k+1) >= 2
                if pRight > 1
                    seq([pRight-1, pRight]) = seq([pRight, pRight-1]);
                end
            else
                if pLeft < L
                    seq([pLeft, pLeft+1]) = seq([pLeft+1, pLeft]);
                end
            end
        end
        pos0 = find(seq == 0);
        segLen = [pos0(1)-1, pos0(2:end)-pos0(1:end-1)-1, L - pos0(end)];
    end
    repaired = seq;
end

function repaired = SegmentBasedCrossover_repairSegments(individual, n, m)
% 辅助函数：修复段内商家-客户对完整性
% 确保每个被0分割的子段内都包含完整的商家-客户对（不分割对），且段内商家在客户前
    idx0 = find(individual == 0);
    if length(idx0) ~= m - 1
        repaired = individual;
        return;
    end
    idx0 = [0, idx0, length(individual)+1];
    segs = cell(m, 1);
    for k = 1:m
        segs{k} = individual(idx0(k)+1 : idx0(k+1)-1);
    end
    
    for i = 1:n
        merchant = i;
        customer = i + n;
        seg_merchant = [];
        seg_customer = [];
        for k = 1:m
            if any(segs{k} == merchant), seg_merchant = k; end
            if any(segs{k} == customer), seg_customer = k; end
        end
        if isempty(seg_merchant), continue; end
        if ~isempty(seg_customer) && seg_customer ~= seg_merchant
            % 将客户点从 seg_customer 移到 seg_merchant（紧跟对应商家后）
            segs{seg_customer} = segs{seg_customer}(segs{seg_customer} ~= customer);
            pos_merchant = find(segs{seg_merchant} == merchant, 1);
            segs{seg_merchant} = [segs{seg_merchant}(1:pos_merchant), customer, segs{seg_merchant}(pos_merchant+1:end)];
        end
    end
    
    repaired = [];
    for k = 1:m
        repaired = [repaired, segs{k}];
        if k < m
            repaired = [repaired, 0];
        end
    end
    % 若首段或末段为空会导致首/尾为 0，统一修正
    repaired = SegmentBasedCrossover_ensureNoZeroAtEnds(repaired);
    % 确保每段至少 2 个基因（至少一对商家-客户）
    repaired = SegmentBasedCrossover_repairSegmentSizes(repaired, m);
end

function safe = SegmentBasedCrossover_isSafeCrossoverPoint(pos_map, k, n)
% 辅助函数：检查交叉点是否安全
    safe = true;
    for s = 1:n
        c = s + n;
        if pos_map(s) < k && pos_map(c) > k
            safe = false;
            return;
        end
    end
end

function priority = SegmentBasedCrossover_calculatePriority(individual, n)
% 辅助函数：计算商家优先级
    priority = zeros(n, 1);
    for i = 1:length(individual)
        gene = individual(i);
        if gene >= 1 && gene <= n
            priority(gene) = i;  % 位置作为优先级
        end
    end
end

function order = SegmentBasedCrossover_extractMerchantOrder(individual, n)
% 辅助函数：提取商家顺序
    order = [];
    for i = 1:length(individual)
        gene = individual(i);
        if gene >= 1 && gene <= n
            order = [order, gene];
        end
    end
    % 去重
    order = unique(order);
end

function seq = SegmentBasedCrossover_ensureNoZeroAtEnds(seq)
% 辅助函数：确保序列首尾不为 0（0 仅作段间分隔，首尾必须是商家或客户点）
    L = length(seq);
    while L >= 1 && seq(1) == 0
        idx = find(seq ~= 0, 1);
        if isempty(idx), break; end
        seq([1, idx]) = seq([idx, 1]);
    end
    while L >= 1 && seq(L) == 0
        idx = find(seq(1:L) ~= 0, 1, 'last');
        if isempty(idx), break; end
        seq([idx, L]) = seq([L, idx]);
    end
end

function repaired = SegmentBasedCrossover_repairSolution(individual, n, m, varDim)
% 辅助函数：修复解的完整性
    repaired = individual;
    
    % 1. 确保包含所有必需的点（1-2*n）
    all_required = 1:(2*n);
    present = unique(repaired(repaired ~= 0));
    missing = setdiff(all_required, present);
    
    if ~isempty(missing)
        % 重复位置：每个重复值的第二次及以后出现的位置（保留第一次）
        duplicate_positions = [];
        for i = 1:length(repaired)
            if repaired(i) ~= 0
                first_pos = find(repaired == repaired(i), 1);
                if first_pos < i
                    duplicate_positions = [duplicate_positions, i];
                end
            end
        end
        % 多余0的位置：只保留前 m-1 个0，其余可用来填缺失
        zero_positions = find(repaired == 0);
        needed_zeros = m - 1;
        extra_zero_positions = zero_positions(needed_zeros+1:end);
        replace_positions = [duplicate_positions, extra_zero_positions];
        
        for i = 1:min(length(missing), length(replace_positions))
            repaired(replace_positions(i)) = missing(i);
        end
    end
    
    % 2. 确保0分隔符数量正确
    zero_count = sum(repaired == 0);
    needed_zeros = m - 1;
    
    if zero_count < needed_zeros
        % 添加缺少的0分隔符，且保证每段至少 2 个基因（不在首尾插 0，且段长>=2）
        non_zero_pos = find(repaired ~= 0);
        N = length(non_zero_pos);
        num_to_add = needed_zeros - zero_count;
        if N >= 2*(num_to_add+1) && num_to_add > 0
            % 插入点选在 non_zero_pos 的索引上：第1段至少2个，中间段至少2个，最后一段至少2个
            % 即 insert_idx 在 [2, N-2] 内且两两间隔至少 2
            first_idx = 2;
            last_idx = N - 2;
            if last_idx >= first_idx
                insert_idx = round(linspace(first_idx, last_idx, num_to_add));
                insert_idx = max(insert_idx, first_idx + (0:num_to_add-1)*2);
                insert_idx = min(insert_idx, last_idx - (num_to_add-1:-1:0)*2);
                for i = length(insert_idx):-1:1
                    pos = non_zero_pos(insert_idx(i));
                    repaired = [repaired(1:pos), 0, repaired(pos+1:end)];
                end
            end
        end
    elseif zero_count > needed_zeros
        % 删除多余的0分隔符
        zero_positions = find(repaired == 0);
        to_remove = zero_positions(needed_zeros+1:end);
        repaired(to_remove) = [];
    end
    
    % 3. 确保维度正确
    if length(repaired) < varDim
        repaired = [repaired, zeros(1, varDim - length(repaired))];
    elseif length(repaired) > varDim
        repaired = repaired(1:varDim);
    end
    
    % 4. 确保首尾不为 0：0 仅作段间分隔，首尾必须是商家或客户点（保证前后都有完整对）
    repaired = SegmentBasedCrossover_ensureNoZeroAtEnds(repaired);
    % 5. 确保每段至少 2 个基因（至少一对商家-客户），避免 "0 18" 这种末段只有一点
    repaired = SegmentBasedCrossover_repairSegmentSizes(repaired, m);
end

function repaired = SegmentBasedCrossover_repairOrder(individual, n)
% 辅助函数：修复顺序约束
    repaired = individual;
    
    for s = 1:n
        c = s + n;
        pos_s = find(repaired == s, 1);
        pos_c = find(repaired == c, 1);
        
        if ~isempty(pos_s) && ~isempty(pos_c) && pos_s > pos_c
            % 将商家移动到客户点之前（删除后客户点位置仍为 pos_c，因 pos_c < pos_s）
            gene = repaired(pos_s);
            repaired(pos_s) = [];
            repaired = [repaired(1:pos_c-1), gene, repaired(pos_c:end)];
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
