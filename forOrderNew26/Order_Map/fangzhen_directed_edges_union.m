function E = fangzhen_directed_edges_union(routes)
% 从多条有向路线（0=配送中心, 1..10 对应 order_data_fangzhen 中 dotss 第 1..10 行）合并有向边并去重
% routes: 元胞数组，每个元素为行向量，形如 [0, u, v, ..., 0]
    if isempty(routes)
        E = zeros(0, 2);
        return;
    end
    E = [];
    for ir = 1:numel(routes)
        seq = routes{ir}(:).';
        if numel(seq) < 2
            continue;
        end
        for k = 1:numel(seq) - 1
            e = [seq(k), seq(k+1)];
            if ~any(ismember(E, e, 'rows'))
                E = [E; e]; %#ok<AGROW>
            end
        end
    end
end
