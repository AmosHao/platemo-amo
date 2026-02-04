function repaired = OrderConstraintRepair(individual, n, m)
% 顺序约束修复函数
% 修复违反顺序约束的解：商家i必须在客户点(i+n)之前访问
%
% 输入:
%   individual - 待修复的个体
%   n - 客户点数量（商家数量）
%   m - 无人机数量
%
% 输出:
%   repaired - 修复后的个体（保持原始维度）

    repaired = individual;
    originalDim = length(individual);
    
    % 找到0分隔符位置，将解分段
    idx0 = find(repaired == 0);
    
    % 确保有 m-1 个0分隔符（如果不足，在末尾补充）
    if length(idx0) < m - 1
        numZerosNeeded = (m - 1) - length(idx0);
        repaired = [repaired, zeros(1, numZerosNeeded)];
        idx0 = find(repaired == 0);
    elseif length(idx0) > m - 1
        % 如果0的数量过多，只保留前m-1个
        idx0 = idx0(1:m-1);
    end
    
    idx0 = [0 idx0 length(repaired)+1];
    
    % 检查每个商家-客户点对
    for merchant_id = 1:n
        customer_id = merchant_id + n;  % 客户点编号
        
        merchant_pos = find(repaired == merchant_id);
        customer_pos = find(repaired == customer_id);
        
        % 如果路径中包含客户点
        if ~isempty(customer_pos)
            % 如果路径中不包含对应商家，跳过（由目标函数惩罚处理）
            if isempty(merchant_pos)
                continue;
            end
            
            % 找到客户点和商家所在的段
            customer_seg = 0;
            merchant_seg = 0;
            for seg_idx = 1:length(idx0)-1
                seg_range = (idx0(seg_idx)+1):(idx0(seg_idx+1)-1);
                if any(repaired(seg_range) == customer_id)
                    customer_seg = seg_idx;
                end
                if any(repaired(seg_range) == merchant_id)
                    merchant_seg = seg_idx;
                end
            end
            
            % 如果商家和客户点在同一段，检查顺序
            if merchant_seg == customer_seg && merchant_seg > 0
                seg_range = (idx0(merchant_seg)+1):(idx0(merchant_seg+1)-1);
                merchant_pos_in_seg = find(repaired(seg_range) == merchant_id, 1);
                customer_pos_in_seg = find(repaired(seg_range) == customer_id, 1);
                
                % 如果商家在客户点之后，交换位置
                if merchant_pos_in_seg >= customer_pos_in_seg
                    segVals = repaired(seg_range);
                    merchant_abs_pos = seg_range(merchant_pos_in_seg);
                    customer_abs_pos = seg_range(customer_pos_in_seg);
                    
                    % 交换位置
                    temp = segVals(merchant_pos_in_seg);
                    segVals(merchant_pos_in_seg) = segVals(customer_pos_in_seg);
                    segVals(customer_pos_in_seg) = temp;
                    
                    repaired(seg_range) = segVals;
                end
            elseif merchant_seg > customer_seg && customer_seg > 0
                % 商家段在客户点段之后，将商家移动到客户点所在段
                seg_range_customer = (idx0(customer_seg)+1):(idx0(customer_seg+1)-1);
                seg_range_merchant = (idx0(merchant_seg)+1):(idx0(merchant_seg+1)-1);
                
                merchant_abs_pos = seg_range_merchant(find(repaired(seg_range_merchant) == merchant_id, 1));
                
                % 移除商家
                repaired(merchant_abs_pos) = [];
                
                % 更新分隔符位置
                idx0_new = find(repaired == 0);
                idx0_new = [0 idx0_new length(repaired)+1];
                
                % 插入到客户点所在段的客户点之前
                seg_range_customer_new = (idx0_new(customer_seg)+1):(idx0_new(customer_seg+1)-1);
                customer_pos_new = seg_range_customer_new(find(repaired(seg_range_customer_new) == customer_id, 1));
                
                repaired = [repaired(1:customer_pos_new-1) merchant_id repaired(customer_pos_new:end)];
            end
        end
    end
    
    % 确保维度一致：如果长度改变，调整到原始维度
    if length(repaired) ~= originalDim
        if length(repaired) < originalDim
            % 如果太短，在末尾补0
            repaired = [repaired, zeros(1, originalDim - length(repaired))];
        else
            % 如果太长，截断
            repaired = repaired(1:originalDim);
        end
    end
end
