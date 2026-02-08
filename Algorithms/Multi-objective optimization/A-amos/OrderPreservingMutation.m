function offspring = OrderPreservingMutation(individual, m, varDim, n)
% 保持顺序约束的变异操作
% 实现三种变异方式：
% 1. 把某个无人机的飞行段中的某个商家移到最后
% 2. 某段的某个客户点移到首位
% 3. 移动0分隔符的位置（必须包含完整的商家-客户点对）
%
% 输入:
%   individual - 待变异的个体
%   m - 无人机数量（段数）
%   varDim - 变量维度
%   n - 商家数量（用于确定必需的点：1到2*n）
%
% 输出:
%   offspring - 变异后的个体（保证：m-1个0分隔符，1-2*n每个值只出现一次，满足顺序约束）

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
    end
    
    % 分段：根据0分隔符将路径分成m段
    idx0_extended = [0, idx0, length(individual)+1];
    segments = cell(m, 1);
    for k = 1:m
        segments{k} = individual(idx0_extended(k)+1 : idx0_extended(k+1)-1);
    end
    
    % 随机选择一种变异方式
    mutation_type = randi(3);
    
    switch mutation_type
        case 1
            % 方式1：把某个无人机的飞行段中的某个商家移到最后
            offspring = mutation_type1(individual, segments, idx0_extended, m, n);
            
        case 2
            % 方式2：某段的某个客户点移到首位
            offspring = mutation_type2(individual, segments, idx0_extended, m, n);
            
        case 3
            % 方式3：移动0分隔符的位置（必须包含完整的商家-客户点对）
            offspring = mutation_type3(individual, segments, idx0_extended, m, n, varDim);
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
end

function offspring = mutation_type1(individual, segments, idx0_extended, m, n)
    % 方式1：把某个无人机的飞行段中的某个商家移到最后
    
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
    
    % 将该商家移到该段的最后
    new_seg = seg;
    new_seg(merchant_idx_in_seg) = [];  % 删除原位置
    new_seg = [new_seg, merchant_id];   % 添加到末尾
    
    % 重建解
    offspring = individual;
    start_idx = idx0_extended(seg_idx) + 1;
    end_idx = idx0_extended(seg_idx + 1) - 1;
    offspring(start_idx:end_idx) = new_seg;
end

function offspring = mutation_type2(individual, segments, idx0_extended, m, n)
    % 方式2：某段的某个客户点移到首位
    
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
    
    % 将该客户点移到该段的首位
    new_seg = seg;
    new_seg(customer_idx_in_seg) = [];  % 删除原位置
    new_seg = [customer_id, new_seg];   % 添加到首位
    
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
    % 返回：新的0分隔符位置（在individual中的索引），如果无法移动返回-1
    
    if isempty(seg_before)
        new_zero_pos = -1;
        return;
    end
    
    % seg_before在individual中的位置范围
    seg_start = idx0_extended(zero_idx) + 1;
    seg_end = idx0_extended(zero_idx + 1) - 1;
    
    % 从seg_before末尾往前遍历，找到完整的配对
    % 记录已找到的商家和客户点
    merchants_found = [];  % 记录找到的商家ID
    merchants_pos = [];   % 记录找到的商家在seg_before中的位置
    customers_found = []; % 记录找到的客户点ID
    customers_pos = [];   % 记录找到的客户点在seg_before中的位置
    
    for i = length(seg_before):-1:1
        point_id = seg_before(i);
        
        if point_id >= 1 && point_id <= n
            % 找到商家
            merchant_id = point_id;
            customer_id = merchant_id + n;
            
            % 检查这个商家对应的客户点是否已经找到（在当前位置之前）
            customer_idx = find(customers_found == customer_id);
            if ~isempty(customer_idx)
                % 找到了一个完整的配对（商家在i，客户点在customers_pos(customer_idx)）
                % 计算客户点在individual中的位置
                customer_pos_in_individual = seg_start + customers_pos(customer_idx) - 1;
                
                % 检查是否可以移动：客户点之前不是0或不是头部
                if customer_pos_in_individual > 1 && individual(customer_pos_in_individual - 1) ~= 0
                    % 可以移动到客户点之前
                    new_zero_pos = customer_pos_in_individual - 1;
                    return;
                end
            else
                % 客户点还没找到，记录这个商家
                merchants_found = [merchant_id, merchants_found];
                merchants_pos = [i, merchants_pos];
            end
            
        elseif point_id > n && point_id <= 2*n
            % 找到客户点
            customer_id = point_id;
            merchant_id = customer_id - n;
            
            % 检查这个客户点对应的商家是否已经找到（在当前位置之前）
            merchant_idx = find(merchants_found == merchant_id);
            if ~isempty(merchant_idx)
                % 找到了一个完整的配对（商家在merchants_pos(merchant_idx)，客户点在i）
                % 计算客户点在individual中的位置
                customer_pos_in_individual = seg_start + i - 1;
                
                % 检查是否可以移动：客户点之前不是0或不是头部
                if customer_pos_in_individual > 1 && individual(customer_pos_in_individual - 1) ~= 0
                    % 可以移动到客户点之前
                    new_zero_pos = customer_pos_in_individual - 1;
                    return;
                end
            else
                % 商家还没找到，记录这个客户点
                customers_found = [customer_id, customers_found];
                customers_pos = [i, customers_pos];
            end
        end
    end
    
    % 如果遍历完整个seg_before都没找到可以移动的位置，返回-1
    new_zero_pos = -1;
end

function new_zero_pos = traverse_backward(individual, idx0_extended, zero_idx, seg_after, n)
    % 往后遍历（向seg_after方向）
    % 从0分隔符位置往后遍历seg_after（从开头往后）
    % 逻辑：先遇到客户点，然后往后检查商家的位置
    % 从当前到商家的位置之间遍历过的点，如果没有缺失配对，就可以移动
    % 找到第一个完整配对后，检查商家点之后是否可以移动（不是0且不是尾部）
    % 返回：新的0分隔符位置（在individual中的索引），如果无法移动返回-1
    
    if isempty(seg_after)
        new_zero_pos = -1;
        return;
    end
    
    % seg_after在individual中的位置范围
    seg_start = idx0_extended(zero_idx + 1) + 1;
    seg_end = idx0_extended(zero_idx + 2) - 1;
    
    % 从seg_after开头往后遍历，找到完整的配对
    % 记录已找到的商家和客户点
    merchants_found = [];  % 记录找到的商家ID
    merchants_pos = [];   % 记录找到的商家在seg_after中的位置
    customers_found = []; % 记录找到的客户点ID
    customers_pos = [];   % 记录找到的客户点在seg_after中的位置
    
    for i = 1:length(seg_after)
        point_id = seg_after(i);
        
        if point_id >= 1 && point_id <= n
            % 找到商家
            merchant_id = point_id;
            customer_id = merchant_id + n;
            
            % 检查这个商家对应的客户点是否已经找到（在当前位置之前）
            customer_idx = find(customers_found == customer_id);
            if ~isempty(customer_idx)
                % 找到了一个完整的配对（商家在i，客户点在customers_pos(customer_idx)）
                % 检查从客户点到商家点之间的点是否都是完整的配对
                customer_pos_in_seg = customers_pos(customer_idx);
                segment_between = seg_after(customer_pos_in_seg:i);
                if is_complete_pairs_segment(segment_between, n)
                    % 计算商家点在individual中的位置
                    merchant_pos_in_individual = seg_start + i - 1;
                    
                    % 检查是否可以移动：商家点之后不是0且不是尾部
                    if merchant_pos_in_individual < length(individual) && ...
                       individual(merchant_pos_in_individual + 1) ~= 0
                        % 可以移动到商家点之后
                        new_zero_pos = merchant_pos_in_individual;
                        return;
                    end
                end
            else
                % 客户点还没找到，记录这个商家
                merchants_found = [merchant_id, merchants_found];
                merchants_pos = [i, merchants_pos];
            end
            
        elseif point_id > n && point_id <= 2*n
            % 找到客户点
            customer_id = point_id;
            merchant_id = customer_id - n;
            
            % 检查这个客户点对应的商家是否在后面
            merchant_idx_in_seg = find(seg_after == merchant_id);
            merchant_idx_in_seg = merchant_idx_in_seg(merchant_idx_in_seg > i);
            
            if ~isempty(merchant_idx_in_seg)
                % 找到了一个完整的配对（客户点在i，商家在merchant_idx_in_seg(1)）
                % 计算商家点在individual中的位置
                merchant_pos_in_individual = seg_start + merchant_idx_in_seg(1) - 1;
                
                % 检查从客户点到商家点之间的点是否都是完整的配对
                % 简化：如果客户点和商家点之间没有其他点，或者都是完整的配对，可以移动
                segment_between = seg_after(i:merchant_idx_in_seg(1));
                if is_complete_pairs_segment(segment_between, n)
                    % 检查是否可以移动：商家点之后不是0且不是尾部
                    if merchant_pos_in_individual < length(individual) && ...
                       individual(merchant_pos_in_individual + 1) ~= 0
                        % 可以移动到商家点之后
                        new_zero_pos = merchant_pos_in_individual;
                        return;
                    end
                end
            else
                % 商家还没找到，记录这个客户点
                customers_found = [customer_id, customers_found];
                customers_pos = [i, customers_pos];
            end
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
