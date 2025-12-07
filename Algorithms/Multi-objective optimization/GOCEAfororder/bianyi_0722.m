function child = bianyi_0722(child)
% child 为行向量，例如 [2 1 0 5 3 4 0 6]
% 仅用 0 作为分隔符，首尾无 0

    % 1. 按 0 切分得到若干子路径（不含 0）
    segments = split_by_zero(child);     % cell 数组，每个 cell 是一段
    genes    = [segments{:}];            % 把所有基因拉平成一个序列

    % 2. 在 genes 中随机选取两个不同位置交换
    len = length(genes);
    if len > 1
        pos1 = randi(len);
        pos2 = randi(len);
        while pos2 == pos1
            pos2 = randi(len);
        end
        genes([pos1, pos2]) = genes([pos2, pos1]);
    end

    % 3. 再随机选一段子路径做“插入”变异
    numSeg = numel(segments);
    if numSeg >= 2
        segID = randi(numSeg);          % 随机选一段
        seg   = segments{segID};
        if numel(seg) >= 2
            % 随机选该段内一个基因挪到最前面
            idx   = randi(numel(seg));
            gene  = seg(idx);
            seg   = [gene, seg([1:idx-1, idx+1:end])];
            segments{segID} = seg;
        end
    end

    % 4. 重新用 0 把各段拼接回去
    child = strjoin(cellfun(@num2str, segments, 'UniformOutput', false), ' 0 ');
    child = str2num(child);             % 得到数字行向量
end


% ---------- 子函数：按 0 切分 ----------
function segs = split_by_zero(chrom)
idx = find(chrom == 0);
if isempty(idx)
    segs = {chrom};          % 整条染色体没有 0
    return
end
start = 1;
segs  = cell(1, numel(idx)+1);
for k = 1:numel(idx)
    segs{k} = chrom(start : idx(k)-1);
    start   = idx(k) + 1;
end
segs{end} = chrom(start:end);
segs(cellfun('isempty', segs)) = [];
end