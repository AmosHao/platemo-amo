function R = read_ablation_run_by_obj_rank(xlsxPath, whichObj, varargin)
%READ_ABLATION_RUN_BY_OBJ_RANK 读取消融 run 的 xlsx，按指定目标排序并返回索引与数据
%
% 默认对应 save_single_run_results 写出的表：obj_data, dec_data, clust_data, details_data
%
% 用法:
%   R = read_ablation_run_by_obj_rank;  % 使用默认文件路径、按 f1 升序
%   R = read_ablation_run_by_obj_rank(path, 2);         % 按第 2 个目标(通常为 f2)升序
%   R = read_ablation_run_by_obj_rank(path, 1, 'SortDir', 'descend');
%
% 输入:
%   xlsxPath  - 完整路径；若缺省/空，用下方默认 ablation 文件
%   whichObj  - 1, 2, ...  对应 Obj 的列下标(第 k 个目标)
% 可选 名称-值:
%   'SortDir'     - 'ascend'(默认) 或 'descend'
%   'ReadClust'   - true(默认) 是否读 clust_data
%   'ReadDetails' - true(默认) 是否读 details_data(与 obj 行一一对应)
%   'ReadHV'      - false(默认) 是否把 hv_data/igd_data 一起读进 R
%
% 输出 R 为 struct，主要字段:
%   .nSol            - 解的个数
%   .whichObjective - 排序所用目标列
%   .sortDir        - 'ascend' 或 'descend'
%   .sortPerm       - 列向量，length nSol；sortPerm(k) 为“排序后第 k 个解”在**原始文件行序**中的行号(1..nSol)
%   .rank           - 列向量 1..nSol；1 表示在指定目标下第 1 好(在 SortDir 意义下)
%   .Obj            - 按该目标排好序的 Obj (nSol x M)
%   .Dec            - 同上行的 Dec
%   .ObjUnsorted, .DecUnsorted - 与 Excel 中原始行序一致(未重排)
%   .clustTag        - 若读入；与 .Obj 行对齐(已重排) - 对应该排序下的聚类
%   .details         - 若读入；table，行与 .Obj 对齐(已重排)
%   .file, .sheets  - 路径与已读取的 sheet 列表
%
% 说明: “排名” 指按**单列目标**做 sort；多目标 Pareto 的支配序与此不同。若多目标为最小化
%   问题，一般使用默认 ascend(值越小越靠前)。

    if nargin < 1 || isempty(xlsxPath) || (isstring(xlsxPath) && (strlength(xlsxPath) == 0)) || (ischar(xlsxPath) && numel(xlsxPath) == 0)
        xlsxPath = "D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\ablation\fangzhen2_NSGAIII_G\fangzhen2_NSGAIII_G_run01.xlsx";
    else
        xlsxPath = char(xlsxPath);
    end
    if nargin < 2 || isempty(whichObj)
        whichObj = 1;
    end

    p = inputParser;
    addParameter(p, 'SortDir', 'ascend', @(s) ismember(lower(string(s)), ["ascend", "descend"]));
    addParameter(p, 'ReadClust', true, @islogical);
    addParameter(p, 'ReadDetails', true, @islogical);
    addParameter(p, 'ReadHV', false, @islogical);
    parse(p, varargin{:});
    sortDir = char(p.Results.SortDir);

    if ~isfile(xlsxPath)
        error('read_ablation_run_by_obj_rank:FileNotFound', '找不到文件: %s', xlsxPath);
    end

    objSheet = 'obj_data';
    decSheet = 'dec_data';
    R.file = xlsxPath;
    R.sheets = {objSheet, decSheet};

    Obj = readmatrix(xlsxPath, 'Sheet', objSheet);
    Dec = readmatrix(xlsxPath, 'Sheet', decSheet);
    nSol = size(Obj, 1);
    if size(Dec, 1) ~= nSol
        error('read_ablation_run_by_obj_rank:SizeMismatch', ...
            'obj_data 行数 %d 与 dec_data 行数 %d 不一致', nSol, size(Dec, 1));
    end

    mObj = size(Obj, 2);
    if whichObj < 1 || whichObj > mObj
        error('read_ablation_run_by_obj_rank:BadObjIndex', ...
            'whichObj=%d 超出范围 [1,%d](obj_data 列数)。', whichObj, mObj);
    end

    [~, perm] = sort(Obj(:, whichObj), 1, sortDir);
    perm = perm(:);

    R.nSol = nSol;
    R.whichObjective = whichObj;
    R.sortDir = sortDir;
    R.sortPerm = perm;
    R.rank = (1:nSol)';
    R.Obj = Obj(perm, :);
    R.Dec = Dec(perm, :);
    R.ObjUnsorted = Obj;
    R.DecUnsorted = Dec;

    % .sortPerm(k) 即：排序后第 k 个解 在 原始 obj/dec 行序中的 行号(1..nSol)

    if p.Results.ReadClust
        csh = 'clust_data';
        if isSheetPresent(xlsxPath, csh)
            C = readmatrix(xlsxPath, 'Sheet', csh, 'NumHeaderLines', 1);
            if size(C, 1) == nSol && size(C, 2) >= 3
                R.clustTag = C(perm, 3);
            elseif size(C, 1) == nSol && size(C, 2) >= 1
                R.clustTag = C(perm, end);
            else
                R.clustTag = [];
                warning('read_ablation_run_by_obj_rank:ClustSize', 'clust_data 行/列与预期不符，已跳过。');
            end
            R.sheets{end+1} = csh; %#ok<AGROW>
        else
            R.clustTag = [];
        end
    else
        R.clustTag = [];
    end

    if p.Results.ReadDetails
        dsh = 'details_data';
        if isSheetPresent(xlsxPath, dsh)
            try
                try
                    Td = readtable(xlsxPath, 'Sheet', dsh, 'VariableNamingRule', 'preserve');
                catch
                    Td = readtable(xlsxPath, 'Sheet', dsh);
                end
                if height(Td) == nSol
                    R.details = Td(perm, :);
                else
                    R.details = table();
                    warning('read_ablation_run_by_obj_rank:DetailsSize', 'details_data 行数与 obj_data 不一致，未返回 details。');
                end
            catch ME
                R.details = table();
                warning('read_ablation_run_by_obj_rank:DetailsRead', '读取 details_data 失败: %s', ME.message);
            end
            R.sheets{end+1} = dsh; %#ok<AGROW>
        else
            R.details = table();
        end
    else
        R.details = table();
    end

    if p.Results.ReadHV
        if isSheetPresent(xlsxPath, 'hv_data')
            R.hv = readmatrix(xlsxPath, 'Sheet', 'hv_data');
        end
        if isSheetPresent(xlsxPath, 'igd_data')
            R.igd = readmatrix(xlsxPath, 'Sheet', 'igd_data');
        end
    end
end

function tf = isSheetPresent(xlsPath, sh)
    tf = false;
    try
        s = sheetnames(xlsPath);
    catch
        [~, s] = xlsfinfo(xlsPath);
    end
    if isempty(s)
        return;
    end
    if ischar(s)
        s = {s};
    end
    if isstring(s)
        s = cellstr(s);
    end
    tf = any(strcmpi(s, sh));
end
