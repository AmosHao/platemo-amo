function route = repair_route(route)
%REPAIR_ROUTE  修复一条路径编码，使其满足：
%   - 长度为 12
%   - 恰好 2 个 0，且位置在 2..11
%   - 非 0 元素为 1..10，每个恰好出现一次

    varDim = 12;

    % 1. 裁剪 / 补齐长度（防御式，正常情况下已经是 1x12）
    route = route(:)';                     % 强制行向量
    if numel(route) > varDim
        route = route(1:varDim);
    elseif numel(route) < varDim
        route = [route, ones(1,varDim-numel(route))];
    end

    % 2. 先把明显非法的数剪到 [0,10] 范围（0 代表分隔符）
    route(route < 0)  = 0;
    route(route > 10) = 10;

    % 3. 处理 0 的个数：目标是恰好 2 个 0，位置限定在 2..11
    zero_pos = find(route == 0);

    % 3.1 如果 0 太多：只保留前两个，其余暂时标记为 -1
    if numel(zero_pos) > 2
        keep = zero_pos(1:2);
        drop = zero_pos(3:end);
        route(drop) = -1;
        zero_pos = keep;
    end

    % 3.2 如果 0 不足：在 2..11 范围内补足
    if numel(zero_pos) < 2
        candidates = 2:11;
        candidates = setdiff(candidates, zero_pos);
        need = 2 - numel(zero_pos);
        if numel(candidates) >= need
            add_pos = candidates(randperm(numel(candidates), need));
        else
            add_pos = candidates;
        end
        route(add_pos) = 0;
        zero_pos = [zero_pos, add_pos];
    end

    % 再次限定 0 的位置在 2..11 之间
    zero_pos(zero_pos < 2 | zero_pos > 11) = [];
    if numel(zero_pos) ~= 2
        % 如果仍不满足，就强制重置 0 位
        zero_pos = randperm(10,2) + 1;
    end

    % 4. 修复非 0 部分：必须是 1..10 且不重复
    non_zero_pos = setdiff(1:varDim, zero_pos);
    genes = route(non_zero_pos);

    % 4.1 先清除不在 1..10 范围的，标记为 -1
    genes(genes < 1 | genes > 10) = -1;

    % 4.2 保留每个客户编号第一次出现，多余的设为 -1
    used = false(1,10);
    for i = 1:numel(genes)
        g = genes(i);
        if g >= 1 && g <= 10
            if ~used(g)
                used(g) = true;
            else
                genes(i) = -1;
            end
        end
    end

    % 4.3 把缺失的客户编号填到标记为 -1 的位置
    missing = find(~used);           % 还没出现的客户
    holes   = find(genes == -1);     % 需要填补的位置
    fill_cnt = min(numel(missing), numel(holes));
    if fill_cnt > 0
        genes(holes(1:fill_cnt)) = missing(1:fill_cnt);
    end

    % 如果仍然不够 10 个客户（极端情况），用剩余的随机补齐
    current = genes(genes >= 1 & genes <= 10);
    if numel(current) < 10
        remain = setdiff(1:10, current);
        extra_holes = find(genes == -1);
        fill_cnt = min(numel(remain), numel(extra_holes));
        if fill_cnt > 0
            genes(extra_holes(1:fill_cnt)) = remain(1:fill_cnt);
        end
    end

    % 5. 最终组装：写回 route
    route = zeros(1,varDim);
    route(zero_pos) = 0;
    route(non_zero_pos) = genes;

    % 6. 最后保险：如果仍不合法，完全随机生成一条合法解
    if numel(find(route==0)) ~= 2 || numel(unique(route(route~=0))) ~= 10
        % 随机选择10个客户编号（1-10）
        customers = randperm(10,10);
        % 随机选择2个0的位置（2-11）
        zero_positions = randperm(10,2) + 1;
        route = zeros(1,12);
        route(zero_positions) = 0;
        non_zero_pos = setdiff(1:12, zero_positions);
        route(non_zero_pos) = customers;
    end
end



