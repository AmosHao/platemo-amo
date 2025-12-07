function child = bianyi_0831(parent)
    cust = parent(parent~=0);          % 10 个客户
    zero = find(parent==0);            % 2 个 0 位

    % 确保输入正确
    if length(cust) ~= 10 || length(zero) ~= 2
        % 如果输入不正确，生成一个有效的随机解
        child = generate_valid_solution();
        return;
    end

    % 1) 客户序列变异：增加多种变异方式
    if numel(cust) >= 2 && rand < 0.9  % 提高变异概率
        mutation_type = randi(4);       % 4种变异方式
        
        switch mutation_type
            case 1  % 交换变异
                i = randi(numel(cust)-1); 
                j = randi([i+1 numel(cust)]);
                cust([i j]) = cust([j i]);
            case 2  % 插入变异
                i = randi(numel(cust)); 
                j = randi(numel(cust));
                if i ~= j
                    temp = cust(i);
                    % 删除第i个元素
                    cust = [cust(1:i-1), cust(i+1:end)];
                    % 在第j个位置插入
                    if j <= length(cust)
                        cust = [cust(1:j-1), temp, cust(j:end)];
                    else
                        cust = [cust, temp];
                    end
                end
            case 3  % 倒序变异
                if rand < 0.3
                    start_pos = randi([1, numel(cust)-1]);
                    end_pos = randi([start_pos+1, numel(cust)]);
                    cust(start_pos:end_pos) = fliplr(cust(start_pos:end_pos));
                end
            case 4  % 随机重排
                if rand < 0.2
                    cust = cust(randperm(numel(cust)));
                end
        end
    end

    % 2) 0 位变异：增加移动范围和概率
    if rand < 0.5  % 提高0位变异概率
        for k = 1:2  % 对每个0位都进行变异
            if rand < 0.6  % 每个0位的变异概率
                % 增加移动范围：可以移动到更远的位置
                move_range = randi([-2, 2]);  % 从-1,1改为-2,2
                new_pos = zero(k) + move_range;
                new_pos = max(2, min(11, new_pos));  % 保证在 2..11
                
                % 检查新位置是否与另一个0位冲突
                if new_pos ~= zero(k) && new_pos ~= zero(3-k)
                    zero(k) = new_pos;
                end
            end
        end
    end

    % 3) 重新组装：确保返回1×12的行向量
    child = zeros(1, 12);             % 初始化1×12的零向量
    child(zero) = 0;                   % 在0位置插入0
    % 将客户编号按顺序填入非零位置
    non_zero_pos = setdiff(1:12, zero);
    
    % 确保cust是10个元素，non_zero_pos是10个位置
    if length(cust) ~= 10
        % 如果cust不是10个元素，重新生成
        cust = randperm(10, 10);
    end
    
    if length(non_zero_pos) ~= 10
        % 如果0位不是2个，重新调整
        zero = randperm(10, 2) + 1;
        child = zeros(1, 12);
        child(zero) = 0;
        non_zero_pos = setdiff(1:12, zero);
    end
    
    % 现在确保维度匹配
    child(non_zero_pos) = cust;
end

% 生成有效解的辅助函数
function sol = generate_valid_solution()
    % 生成1×12的有效解
    sol = zeros(1, 12);
    
    % 随机选择10个客户编号（1-10）
    customers = randperm(10, 10);
    
    % 随机选择2个0的位置（2-11）
    zero_positions = randperm(10, 2) + 1;
    
    % 组装解
    sol(zero_positions) = 0;
    non_zero_pos = setdiff(1:12, zero_positions);
    sol(non_zero_pos) = customers;
end