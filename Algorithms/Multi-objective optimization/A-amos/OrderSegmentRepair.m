function sol = OrderSegmentRepair(sol, n, m, varDim)
% 统一段结构修复：确保每段均为完整商家-客户对（偶数长度），可有两个 0 相邻但左右段必须完整
% 供 A_amos 在变异、段交换后调用，避免出现 "0 0 18" 这类仅一点在段内的情况
%
% 输入: sol - 决策向量, n - 商家数, m - 无人机数, varDim - 维度
% 输出: 修复后的解（长度 varDim，m-1 个 0，每段至少 2 个点）

    if isempty(sol) || length(sol) < 2
        return;
    end
    if length(sol) ~= varDim
        if length(sol) < varDim
            sol = [sol, zeros(1, varDim - length(sol))];
        else
            sol = sol(1:varDim);
        end
    end
    sol = EnsureNoZeroAtEnds(sol);
    sol = RepairSegments(sol, n, m);
    sol = RepairSegmentSizes(sol, m);
    % RepairSegmentSizes 为凑“每段至少 2 点”可能与 0 交换，会把已合并的商-客对拆回两段，再并一次
    sol = RepairSegments(sol, n, m);
    sol = RepairSegmentSizes(sol, m);
    if length(sol) ~= varDim
        if length(sol) < varDim
            sol = [sol, zeros(1, varDim - length(sol))];
        else
            sol = sol(1:varDim);
        end
    end
end

function seq = EnsureNoZeroAtEnds(seq)
    L = length(seq);
    while L >= 1 && seq(1) == 0
        idx = find(seq ~= 0, 1);
        if isempty(idx), break; end
        seq([1, idx]) = seq([idx, 1]);
    end
    while L >= 1 && seq(L) == 0
        if L >= 2 && seq(L-1) == 0
            seqShort = seq(1:L-1);
            nonZeroIdx = find(seqShort ~= 0);
            inserted = false;
            for k = 1:floor((length(nonZeroIdx)-1)/2)
                if 2*k+1 > length(nonZeroIdx), break; end
                afterPos = nonZeroIdx(2*k);
                if afterPos < length(seqShort) && nonZeroIdx(2*k+1) == afterPos + 1
                    seq = [seqShort(1:afterPos), 0, seqShort(afterPos+1:end)];
                    inserted = true;
                    break;
                end
            end
            if ~inserted && length(nonZeroIdx) >= 2
                afterPos = nonZeroIdx(min(2, length(nonZeroIdx)-1));
                seq = [seqShort(1:afterPos), 0, seqShort(afterPos+1:end)];
            else
                break;
            end
        else
            idx = find(seq(1:L) ~= 0, 1, 'last');
            if isempty(idx), break; end
            seq([idx, L]) = seq([L, idx]);
        end
        L = length(seq);
    end
end

function repaired = RepairSegmentSizes(seq, m)
    pos0 = find(seq == 0);
    if isempty(pos0) || length(pos0) ~= m - 1
        repaired = seq;
        return;
    end
    L = length(seq);
    segLen = [pos0(1)-1, pos0(2:end)-pos0(1:end-1)-1, L - pos0(end)];
    maxIter = L;
    iter = 0;
    while any(segLen < 2) && iter < maxIter
        iter = iter + 1;
        consecutiveZeros = find(diff(pos0) == 1);
        if ~isempty(consecutiveZeros)
            removePos = pos0(consecutiveZeros(1) + 1);
            seqShort = [seq(1:removePos-1), seq(removePos+1:end)];
            nonZeroIdx = find(seqShort ~= 0);
            inserted = false;
            for k = 1:floor((length(nonZeroIdx)-1)/2)
                if 2*k+1 > length(nonZeroIdx), break; end
                afterPos = nonZeroIdx(2*k);
                if afterPos < length(seqShort) && nonZeroIdx(2*k+1) == afterPos + 1
                    seq = [seqShort(1:afterPos), 0, seqShort(afterPos+1:end)];
                    inserted = true;
                    break;
                end
            end
            if ~inserted && length(nonZeroIdx) >= 3
                afterPos = nonZeroIdx(2);
                seq = [seqShort(1:afterPos), 0, seqShort(afterPos+1:end)];
            else
                seq = seqShort;
                if length(nonZeroIdx) >= 2
                    afterPos = nonZeroIdx(min(2, length(nonZeroIdx)-1));
                    seq = [seq(1:afterPos), 0, seq(afterPos+1:end)];
                end
            end
            if length(seq) > L
                seq = seq(1:L);
            end
            pos0 = find(seq == 0);
            segLen = [pos0(1)-1, pos0(2:end)-pos0(1:end-1)-1, L - pos0(end)];
            continue;
        end
        if segLen(1) < 2
            p = pos0(1);
            if p < L && seq(p+1) ~= 0
                seq([p, p+1]) = seq([p+1, p]);
            end
        elseif segLen(m) < 2
            p = pos0(end);
            if p > 1 && seq(p-1) ~= 0
                seq([p-1, p]) = seq([p, p-1]);
            end
        else
            k = find(segLen < 2, 1);
            pLeft = pos0(k-1);
            pRight = pos0(k);
            if k > 1 && segLen(k-1) >= 2 && pLeft < L && seq(pLeft+1) ~= 0
                seq([pLeft, pLeft+1]) = seq([pLeft+1, pLeft]);
            elseif k < m && segLen(k+1) >= 2 && pRight > 1 && seq(pRight-1) ~= 0
                seq([pRight-1, pRight]) = seq([pRight, pRight-1]);
            elseif pLeft < L && seq(pLeft+1) ~= 0
                seq([pLeft, pLeft+1]) = seq([pLeft+1, pLeft]);
            end
        end
        pos0 = find(seq == 0);
        segLen = [pos0(1)-1, pos0(2:end)-pos0(1:end-1)-1, L - pos0(end)];
    end
    repaired = seq;
end

function repaired = RepairSegments(individual, n, m)
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
    repaired = EnsureNoZeroAtEnds(repaired);
    repaired = RepairSegmentSizes(repaired, m);
end
