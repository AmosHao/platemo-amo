function child = genetic_crossover_0831(parent1, parent2)
% parent1/2 : 1×12 向量，含 10 个 1…10 的客户编号 + 2 个 0（2…11 内）

    %% 1. 提取 10 个客户编号
    cust1 = parent1(parent1 ~= 0);   % 10×1
    cust2 = parent2(parent2 ~= 0);   % 10×1

    %% 2. 交叉：增加多种交叉方式
    n = 10;  % 改为10个客户
    crossover_type = randi(3);  % 3种交叉方式
    
    switch crossover_type
        case 1  % 顺序交叉 (OX)
            cp = randi([1 n-1]);
            segLeft = cust1(1:cp);          % 父 1 前段
            segRight = cust2(cp+1:end);    % 父 2 后段
            
            % 去重：把父 1 中已出现的基因从 segRight 中剔除
            segRight = setdiff(segRight, segLeft);
            % 把缺失的客户按顺序补到 segRight 后面，凑够 10 位
            need = setdiff(1:10, [segLeft, segRight]);
            child_cust = [segLeft, segRight, need];
            
        case 2  % 均匀交叉
            mask = rand(1, n) < 0.5;  % 随机掩码
            child_cust = zeros(1, n);
            child_cust(mask) = cust1(mask);
            child_cust(~mask) = cust2(~mask);
            
            % 处理重复值：确保最终有10个不同的客户编号
            while length(unique(child_cust)) < 10
                missing = setdiff(1:10, child_cust);
                if ~isempty(missing)
                    % 随机选择一个重复位置替换
                    duplicates = find(histcounts(child_cust, 1:11) > 1);
                    if ~isempty(duplicates)
                        dup_pos = find(child_cust == duplicates(1), 1);
                        child_cust(dup_pos) = missing(1);
                    end
                end
            end
            
        case 3  % 两点交叉
            cp1 = randi([1, n-2]);
            cp2 = randi([cp1+1, n-1]);
            
            child_cust = cust1;
            child_cust(cp1+1:cp2) = cust2(cp1+1:cp2);
            
            % 处理重复值：确保最终有10个不同的客户编号
            while length(unique(child_cust)) < 10
                missing = setdiff(1:10, child_cust);
                if ~isempty(missing)
                    % 随机选择一个重复位置替换
                    duplicates = find(histcounts(child_cust, 1:11) > 1);
                    if ~isempty(duplicates)
                        dup_pos = find(child_cust == duplicates(1), 1);
                        child_cust(dup_pos) = missing(1);
                    end
                end
            end
    end

    %% 3. 0 位置：增加随机性
    zeroMask = parent1 == 0;           % 12×1 logical
    if rand < 0.7  % 70%概率继承父1的0位
        zeroMask = parent1 == 0;
    elseif rand < 0.5  % 15%概率继承父2的0位
        zeroMask = parent2 == 0;
    else  % 15%概率随机生成新的0位
        zeroMask = false(1, 12);
        zero_pos = randperm(10, 2) + 1;  % 在2-11位置随机选择2个
        zeroMask(zero_pos) = true;
    end

    %% 4. 组装：把 0 插回去，确保返回1×12的行向量
    child = zeros(1, 12);              % 初始化1×12的零向量
    child(zeroMask) = 0;               % 把 2 个 0 塞回对应位
    
    % 确保child_cust是10个元素
    if length(child_cust) ~= 10
        % 如果child_cust不是10个元素，重新生成
        child_cust = randperm(10, 10);
    end
    
    % 将客户编号按顺序填入非零位置
    non_zero_pos = find(~zeroMask);
    
    % 确保non_zero_pos是10个位置
    if length(non_zero_pos) ~= 10
        % 如果0位不是2个，重新调整
        zeroMask = false(1, 12);
        zero_pos = randperm(10, 2) + 1;  % 在2-11位置随机选择2个
        zeroMask(zero_pos) = true;
        child(zeroMask) = 0;
        non_zero_pos = find(~zeroMask);
    end
    
    % 最终检查：确保所有变量都有正确的维度
    if length(child_cust) ~= 10 || length(non_zero_pos) ~= 10
        % 如果仍有问题，生成一个完全随机的解
        child = zeros(1, 12);
        zero_pos = randperm(10, 2) + 1;
        child(zero_pos) = 0;
        non_zero_pos = setdiff(1:12, zero_pos);
        child_cust = randperm(10, 10);
    end
    
    % 现在确保维度匹配
    child(non_zero_pos) = child_cust;
    
    % 最终验证
    if length(child) ~= 12
        error('Generated child has wrong length: %d', length(child));
    end
end