function [offspring, dbgInfo] = OrderPreservingMutation(individual, m, varDim, n, explorationRate)
% 保持顺序约束的变异操作
% 实现四种变异方式：
% 1. 把某个无人机的飞行段中的某个商家移到最后
% 2. 某段的某个客户点移到首位
% 3. 移动0分隔符的位置（必须包含完整的商家-客户点对）
% 4. 大步长：将完整商家-客户对从源段迁移到目标段（跨段移动）
%
% 输入:
%   individual     - 待变异的个体
%   m              - 无人机数量（段数）
%   varDim         - 变量维度
%   n              - 商家数量（用于确定必需的点：1到2*n）
%   explorationRate - T4大步长变异的概率（默认0.25，停滞时传入更高值）
%
% 输出:
%   offspring - 变异后的个体（保证：m-1个0分隔符，1-2*n每个值只出现一次，满足顺序约束）
%   dbgInfo   - 调试信息结构体（strategy, cloneDefended）

    if nargin < 5 || isempty(explorationRate)
        explorationRate = 0.25;
    end
    dbgInfo.strategy      = 0;
    dbgInfo.cloneDefended = false;

    % 确保输入维度正确
    if length(individual) ~= varDim
        fprintf('Warning: 维度不正确');
        return;
    end
    
    % 找到0分隔符的位置
    idx0 = find(individual == 0);
    
    % 如果0分隔符数量不对，先修正
    if length(idx0) ~= m - 1
        fprintf('Warning: 0分隔符数量错误');
        return;
    end
    
    % 分段：根据0分隔符将路径分成m段
    idx0_extended = [0, idx0, length(individual)+1];
    segments = cell(m, 1);
    for k = 1:m
        segments{k} = individual(idx0_extended(k)+1 : idx0_extended(k+1)-1);
    end
    
    % 选择变异方式：以 explorationRate 概率执行大步长T4，其余均匀选策略1-3
    % 停滞时外部传入更高的 explorationRate，提升大步长比例
    if rand < explorationRate
        mutation_type = 4;
    else
        mutation_type = randi(3);
    end
    
    dbgInfo.strategy = mutation_type;
    switch mutation_type
        case 1
            offspring = mutation_type1(individual, segments, idx0_extended, m, n);
        case 2
            offspring = mutation_type2(individual, segments, idx0_extended, m, n);
        case 3
            offspring = mutation_type3(individual, segments, idx0_extended, m, n, varDim);
        case 4
            offspring = mutation_type4_segSwap(individual, segments, m, n, varDim);
    end
    
    % 确保维度正确
    if length(offspring) ~= varDim
        if length(offspring) < varDim
            offspring = [offspring, zeros(1, varDim - length(offspring))];
        else
            offspring = offspring(1:varDim);
        end
    end
    
    % 最终验证：确保0分隔符数量正确
    zeroCount = sum(offspring == 0);
    if zeroCount ~= m - 1
        % 如果0分隔符数量不对，尝试修正
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
            else
                offspring = [offspring, zeros(1, (m-1) - zeroCount)];
            end
        elseif zeroCount > m - 1
            zeroIdx = find(offspring == 0);
            offspring(zeroIdx(m:end)) = [];
        end
        % 重新调整维度
        if length(offspring) ~= varDim
            if length(offspring) < varDim
                offspring = [offspring, zeros(1, varDim - length(offspring))];
            else
                offspring = offspring(1:varDim);
            end
        end
    end
    
    % 若变异退化为克隆（子代与父代完全相同），升级为 type4 大步长扰动
    if isequal(offspring, individual)
        dbgInfo.cloneDefended = true;
        offspring = mutation_type4_segSwap(individual, segments, m, n, varDim);
        % 若 type4 也失败（单段问题等），做简单的段内随机对交换
        if isequal(offspring, individual)
            idx0 = find(individual == 0);
            idx0_ext2 = [0, idx0, varDim+1];
            for k = 1:m
                seg = individual(idx0_ext2(k)+1 : idx0_ext2(k+1)-1);
                merchants = seg(seg >= 1 & seg <= n);
                if length(merchants) >= 2
                    sw = randperm(length(merchants), 2);
                    m1 = merchants(sw(1)); c1 = m1 + n;
                    m2 = merchants(sw(2)); c2 = m2 + n;
                    if any(seg == c1) && any(seg == c2)
                        seg_stripped = seg(seg ~= m1 & seg ~= c1 & seg ~= m2 & seg ~= c2);
                        new_seg = [seg_stripped, m1, c1, m2, c2];
                        rebuilt = individual;
                        rebuilt(idx0_ext2(k)+1 : idx0_ext2(k+1)-1) = new_seg(1:length(seg));
                        offspring = rebuilt;
                        break;
                    end
                end
            end
        end
    end
end

function offspring = mutation_type1(individual, segments, idx0_extended, m, n)
    % 方式1：把某个无人机的飞行段中的某个商家移到首位
    
    % 找到包含商家的段
    valid_segments = [];
    for k = 1:m
        seg = segments{k};
        if ~isempty(seg) && any(seg >= 1 & seg <= n)
            valid_segments = [valid_segments, k];
        end
    end
    
    if isempty(valid_segments)
        offspring = individual;
        return;
    end
    
    % 随机选择一个段
    seg_idx = valid_segments(randi(length(valid_segments)));
    seg = segments{seg_idx};
    
    % 找到该段中的商家位置
    merchant_positions = find(seg >= 1 & seg <= n);
    if isempty(merchant_positions)
        offspring = individual;
        return;
    end
    
    % 随机选择一个商家
    merchant_idx_in_seg = merchant_positions(randi(length(merchant_positions)));
    merchant_id = seg(merchant_idx_in_seg);
    
    % 将该商家移到该段的首位
    new_seg = seg;
    new_seg(merchant_idx_in_seg) = [];  % 删除原位置
    new_seg = [merchant_id, new_seg];   % 添加到首位
    
    % 重建解
    offspring = individual;
    start_idx = idx0_extended(seg_idx) + 1;
    end_idx = idx0_extended(seg_idx + 1) - 1;
    offspring(start_idx:end_idx) = new_seg;
end

function offspring = mutation_type2(individual, segments, idx0_extended, m, n)
    % 方式2：某段的某个客户点移到最后
    
    % 找到包含客户点的段
    valid_segments = [];
    for k = 1:m
        seg = segments{k};
        if ~isempty(seg) && any(seg > n & seg <= 2*n)
            valid_segments = [valid_segments, k];
        end
    end
    
    if isempty(valid_segments)
        offspring = individual;
        return;
    end
    
    % 随机选择一个段
    seg_idx = valid_segments(randi(length(valid_segments)));
    seg = segments{seg_idx};
    
    % 找到该段中的客户点位置
    customer_positions = find(seg > n & seg <= 2*n);
    if isempty(customer_positions)
        offspring = individual;
        return;
    end
    
    % 随机选择一个客户点
    customer_idx_in_seg = customer_positions(randi(length(customer_positions)));
    customer_id = seg(customer_idx_in_seg);
    
    % 将该客户点移到该段的最后
    new_seg = seg;
    new_seg(customer_idx_in_seg) = [];  % 删除原位置
    new_seg = [new_seg, customer_id];   % 添加到最后
    
    % 重建解
    offspring = individual;
    start_idx = idx0_extended(seg_idx) + 1;
    end_idx = idx0_extended(seg_idx + 1) - 1;
    offspring(start_idx:end_idx) = new_seg;
end

function offspring = mutation_type3(individual, segments, idx0_extended, m, n, varDim)
    % 方式3：移动0分隔符的位置（必须包含完整的商家-客户点对）
    % 通过遍历找到完整的配对，然后移动0分隔符
    
    if m <= 1
        offspring = individual;
        return;
    end
    
    % 找到可以移动的0分隔符（前后都有元素）
    movable_zeros = [];
    for z = 1:(m-1)
        seg_before = segments{z};
        seg_after = segments{z+1};
        if ~isempty(seg_before) || ~isempty(seg_after)
            movable_zeros = [movable_zeros, z];
        end
    end
    
    if isempty(movable_zeros)
        offspring = individual;
        return;
    end
    
    % 随机选择一个0分隔符
    zero_idx = movable_zeros(randi(length(movable_zeros)));
    
    % 获取该0分隔符前后的段
    seg_before = segments{zero_idx};
    seg_after = segments{zero_idx + 1};
    
    % 决定移动方向（向前或向后）
    direction = randi(2) - 1;  % 0=向前遍历（向seg_before方向），1=向后遍历（向seg_after方向）
    
    new_zero_pos = -1;  % 新的0分隔符位置（在individual中的索引）
    
    if direction == 0 && ~isempty(seg_before)
        % 往前遍历（向seg_before方向）
        % 从0分隔符位置往前遍历seg_before（从末尾往前）
        new_zero_pos = traverse_forward(individual, idx0_extended, zero_idx, seg_before, n);
        
    elseif direction == 1 && ~isempty(seg_after)
        % 往后遍历（向seg_after方向）
        % 从0分隔符位置往后遍历seg_after（从开头往后）
        new_zero_pos = traverse_backward(individual, idx0_extended, zero_idx, seg_after, n);
    end
    
    if new_zero_pos < 0
        % 无法移动
        offspring = individual;
        return;
    end
    
    % 移动0分隔符到新位置
    % 先删除原来的0
    old_zero_pos = idx0_extended(zero_idx + 1);  % 0分隔符在individual中的位置
    temp_individual = individual;
    temp_individual(old_zero_pos) = [];  % 删除原来的0
    
    % 调整新位置（因为删除了一个元素，位置需要调整）
    if new_zero_pos > old_zero_pos
        new_zero_pos = new_zero_pos - 1;
    end
    
    % 在新位置插入0
    if new_zero_pos == 0
        temp_individual = [0, temp_individual];
    elseif new_zero_pos >= length(temp_individual)
        temp_individual = [temp_individual, 0];
    else
        temp_individual = [temp_individual(1:new_zero_pos), 0, temp_individual(new_zero_pos+1:end)];
    end
    
    offspring = temp_individual;
    
    % 确保维度正确
    if length(offspring) ~= varDim
        if length(offspring) < varDim
            offspring = [offspring, zeros(1, varDim - length(offspring))];
        else
            offspring = offspring(1:varDim);
        end
    end
end

function new_zero_pos = traverse_forward(individual, idx0_extended, zero_idx, seg_before, n)
    % 往前遍历（向seg_before方向）
    % 从0分隔符位置往前遍历seg_before（从末尾往前）
    % 只寻找商家，检查已经遍历过的点是否只包含当前已遍历商家对应的客户点
    % 返回：新的0分隔符位置（在individual中的索引），如果无法移动返回-1
    
    if isempty(seg_before)
        new_zero_pos = -1;
        return;
    end
    
    % seg_before在individual中的位置范围
    seg_start = idx0_extended(zero_idx) + 1;
    seg_end = idx0_extended(zero_idx + 1) - 1;
    
    % 从seg_before末尾往前遍历，只关注商家
    % 记录已遍历过的商家和客户点
    merchants_traversed = [];  % 记录已遍历过的商家ID
    customers_traversed = []; % 记录已遍历过的客户点ID
    
    for i = length(seg_before):-1:1
        point_id = seg_before(i);
        
        if point_id >= 1 && point_id <= n
            % 找到商家
            merchant_id = point_id;
            customer_id = merchant_id + n;
            
            % 检查已经遍历过的点（从末尾到当前位置之后）是否只包含【当前+已遍历】商家对应的客户点
            % 已遍历过的点：从i+1到length(seg_before)之间的点
            if i < length(seg_before)
                traversed_points = seg_before(i+1:end);
                traversed_customers = traversed_points(traversed_points > n & traversed_points <= 2*n);
                
                % 当前已遍历过的商家 = 当前商家 + 之前遍历到的商家
                expected_customers = [merchants_traversed, merchant_id] + n;
                
                % 检查已遍历的客户点是否只包含这些商家对应的客户点
                if isempty(traversed_customers) || all(ismember(traversed_customers, expected_customers))
                    % 满足条件：可以在这个商家之前插入0
                    % 计算商家在individual中的位置
                    merchant_pos_in_individual = seg_start + i - 1;
                    
                    % 检查是否可以移动：商家之前不是0且不是头部
                    if merchant_pos_in_individual > 1 && individual(merchant_pos_in_individual - 1) ~= 0
                        % 可以移动到商家之前
                        new_zero_pos = merchant_pos_in_individual - 1;
                        return;
                    end
                end
            end
            
            % 记录这个商家
            merchants_traversed = [merchants_traversed, merchant_id];
            
        elseif point_id > n && point_id <= 2*n
            % 找到客户点，只记录，不处理
            customer_id = point_id;
            customers_traversed = [customers_traversed, customer_id];
        end
    end
    
    % 如果遍历完整个seg_before都没找到可以移动的位置，返回-1
    new_zero_pos = -1;
end

function new_zero_pos = traverse_backward(individual, idx0_extended, zero_idx, seg_after, n)
    % 往后遍历（向seg_after方向）
    % 从0分隔符位置往后遍历seg_after（从开头往后）
    % 只寻找客户点，检查已经遍历过的点是否只包含当前已遍历客户点对应的商家
    % 返回：新的0分隔符位置（在individual中的索引），如果无法移动返回-1
    
    if isempty(seg_after)
        new_zero_pos = -1;
        return;
    end
    
    % seg_after在individual中的位置范围
    seg_start = idx0_extended(zero_idx + 1) + 1;
    seg_end = idx0_extended(zero_idx + 2) - 1;
    
    % 从seg_after开头往后遍历，只关注客户点
    % 记录已遍历过的商家和客户点
    merchants_traversed = [];  % 记录已遍历过的商家ID
    customers_traversed = []; % 记录已遍历过的客户点ID
    
    for i = 1:length(seg_after)
        point_id = seg_after(i);
        
        if point_id > n && point_id <= 2*n
            % 找到客户点
            customer_id = point_id;
            merchant_id = customer_id - n;
            
            % 检查已经遍历过的点（从开头到当前位置之前）是否只包含【当前+已遍历】客户点对应的商家
            % 已遍历过的点：从1到i-1之间的点
            if i > 1
                traversed_points = seg_after(1:i-1);
                traversed_merchants = traversed_points(traversed_points >= 1 & traversed_points <= n);
                
                % 当前已遍历过的客户点 = 当前客户点 + 之前遍历到的客户点
                expected_merchants = [customers_traversed, customer_id] - n;
                
                % 检查已遍历的商家是否只包含这些客户点对应的商家
                if isempty(traversed_merchants) || all(ismember(traversed_merchants, expected_merchants))
                    % 满足条件：可以在这个客户点之后插入0
                    % 计算客户点在individual中的位置
                    customer_pos_in_individual = seg_start + i - 1;
                    
                    % 检查是否可以移动：客户点之后不是0且不是尾部
                    if customer_pos_in_individual < length(individual) && ...
                       individual(customer_pos_in_individual + 1) ~= 0
                        % 可以移动到客户点之后
                        new_zero_pos = customer_pos_in_individual;
                        return;
                    end
                end
            end
            
            % 记录这个客户点
            customers_traversed = [customers_traversed, customer_id];
            
        elseif point_id >= 1 && point_id <= n
            % 找到商家，只记录，不处理
            merchant_id = point_id;
            merchants_traversed = [merchants_traversed, merchant_id];
        end
    end
    
    % 如果遍历完整个seg_after都没找到可以移动的位置，返回-1
    new_zero_pos = -1;
end

function is_complete = is_complete_pairs_segment(segment, n)
    % 检查segment是否只包含完整的商家-客户点对
    % 简化检查：如果segment长度为偶数，且包含的商家和客户点数量匹配，认为是完整的
    
    if isempty(segment)
        is_complete = true;
        return;
    end
    
    merchants = segment(segment >= 1 & segment <= n);
    customers = segment(segment > n & segment <= 2*n);
    
    % 检查每个商家是否都有对应的客户点
    for i = 1:length(merchants)
        merchant_id = merchants(i);
        customer_id = merchant_id + n;
        if ~any(customers == customer_id)
            is_complete = false;
            return;
        end
    end
    
    % 检查每个客户点是否都有对应的商家
    for i = 1:length(customers)
        customer_id = customers(i);
        merchant_id = customer_id - n;
        if ~any(merchants == merchant_id)
            is_complete = false;
            return;
        end
    end
    
    is_complete = true;
end

function offspring = mutation_type4_segSwap(individual, segments, m, n, varDim)
% 大步长变异：从源段随机抽一个完整商家-客户对，移到目标段的随机位置
% 这与 mutation_type3（移动分隔符）的区别在于：
%   - mutation_type3 只移动 0，整对留在原段边界附近
%   - mutation_type4 直接把整对从源段"剪切"到目标段任意插入点，步长更大

    offspring = individual;  % 默认不变
    
    % 找有多于1对的源段（至少含2个商家-客户对）
    srcCandidates = [];
    for k = 1:m
        seg = segments{k};
        merchants = seg(seg >= 1 & seg <= n);
        if length(merchants) >= 2
            srcCandidates(end+1) = k;
        end
    end
    if isempty(srcCandidates), return; end
    
    % 随机选源段
    srcIdx = srcCandidates(randi(length(srcCandidates)));
    srcSeg = segments{srcIdx};
    srcMerchants = srcSeg(srcSeg >= 1 & srcSeg <= n);
    
    % 从源段随机选一个商家-客户对
    chosenMerchant = srcMerchants(randi(length(srcMerchants)));
    chosenCustomer = chosenMerchant + n;
    
    % 确保客户点也在源段中（顺序约束）
    if ~any(srcSeg == chosenCustomer), return; end
    
    % 从源段移除商家和客户点，保持剩余元素相对顺序
    newSrcSeg = srcSeg(srcSeg ~= chosenMerchant & srcSeg ~= chosenCustomer);
    
    % 随机选目标段（不同于源段）
    dstOptions = setdiff(1:m, srcIdx);
    if isempty(dstOptions), return; end
    dstIdx = dstOptions(randi(length(dstOptions)));
    dstSeg = segments{dstIdx};
    
    % 在目标段随机位置插入商家-客户对（商家在前，客户紧随）
    % 插入点范围：0 表示插到最前，length(dstSeg) 表示插到最后
    insertPos = randi(length(dstSeg) + 1) - 1;
    newDstSeg = [dstSeg(1:insertPos), chosenMerchant, chosenCustomer, dstSeg(insertPos+1:end)];
    
    % 更新 segments，重新拼接
    segments{srcIdx} = newSrcSeg;
    segments{dstIdx} = newDstSeg;
    
    % 重建个体：段之间插入 0 分隔符
    rebuilt = [];
    for k = 1:m
        rebuilt = [rebuilt, segments{k}];
        if k < m
            rebuilt = [rebuilt, 0];
        end
    end
    
    % 验证维度和内容
    if length(rebuilt) ~= varDim, return; end
    allRequired = 1:(2*n);
    if ~isequal(sort(rebuilt(rebuilt ~= 0)), sort(allRequired)), return; end
    
    % 验证顺序约束（每个商家在其客户点之前）
    for k = 1:n
        posM = find(rebuilt == k, 1);
        posC = find(rebuilt == k+n, 1);
        if isempty(posM) || isempty(posC) || posM > posC, return; end
    end
    
    offspring = rebuilt;
end
