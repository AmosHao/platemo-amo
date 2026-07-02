%% 读 ablation 单次 run 的 xlsx，支持「全表重排」+「按第几列目标的第几小」取单条
% 表结构需与 save_single_run_results 一致：obj_data, dec_data,（可选）clust_data, details_data
%
% 在 MATLAB 中：addpath( ... 'Order_ComputeFunction' ) 后  运行: read_ablation_run_by_obj_script
%
% 「第 k 小」= 在指定**列**上按**升序**编号：k=1 最小，k=2 次小，…（有并列值时仍给满 nSol 个名次）

%% ========== 参数区（只改这里）==========
xlsxPath   = "D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\ablation\fangzhen7_NSGAIII_G\fangzhen7_NSGAIII_G_run01.xlsx";

% --- 要查的「第 whichCol 列 第 kSmallest 小」的解（原始行索引 + 目标 + dec）---
whichCol   = 1;              % 第几列目标，如 1=f1(能耗J)，2=f2(时间s)
kSmallest  = 1;              % 第几小：1=该列中最小，2=次小，…

% --- 结构体 R：整表按 whichCol+sortDir 重排（与上面「第几小」可独立理解）---
whichObj   = 1;              % 全表排序用的列；可与 whichCol 不同
sortDir    = "ascend";      % "ascend" 或 "descend"

readClust  = true;            % 是否读 clust_data
readDetails= true;            % 是否读 details_data
readHV     = false;            % 是否再读 hv_data / igd_data

%% ========== 以下一般不用动 ==========
xlsxPath   = char(xlsxPath);
sortDir    = char(lower(strtrim(sortDir)));
if ~ismember(sortDir, {'ascend', 'descend'})
    error('sortDir 须为 ascend 或 descend。');
end
if ~isfile(xlsxPath)
    error('找不到文件: %s', xlsxPath);
end

try
    allSheets = sheetnames(xlsxPath);
catch
    [~, allSheets] = xlsfinfo(xlsxPath);
end
if isempty(allSheets)
    error('无法读取 xlsx 的 sheet 名: %s', xlsxPath);
end
if ischar(allSheets)
    allSheets = {allSheets};
end
if isstring(allSheets)
    allSheets = cellstr(allSheets);
end
hasSheet = @(sh) any(strcmpi(allSheets, sh));

objSheet = 'obj_data';
decSheet = 'dec_data';
if ~hasSheet(objSheet) || ~hasSheet(decSheet)
    error('缺少 sheet: %s 或 %s', objSheet, decSheet);
end

Obj = readmatrix(xlsxPath, 'Sheet', objSheet);
Dec = readmatrix(xlsxPath, 'Sheet', decSheet);
nSol  = size(Obj, 1);
nObjC = size(Obj, 2);
if size(Dec, 1) ~= nSol
    error('obj_data 行数 %d 与 dec_data 行数 %d 不一致', nSol, size(Dec, 1));
end
if whichObj < 1 || whichObj > nObjC
    error('whichObj=%d 超出 [1, %d]', whichObj, nObjC);
end
if whichCol < 1 || whichCol > nObjC
    error('whichCol=%d 超出 [1, %d]', whichCol, nObjC);
end
if kSmallest < 1 || kSmallest > nSol || kSmallest ~= floor(kSmallest)
    error('kSmallest 须为 1..%d 的整数。', nSol);
end

% 第 whichCol 列 升序 中的第 k 个 → 在原始表中的行号
[~, permSmall] = sort(Obj(:, whichCol), 1, 'ascend');
idxInSheet     = permSmall(kSmallest);   % 即 obj_data / dec_data 的「行索引」
objK           = Obj(idxInSheet, :);     % 该行全部目标（1×M 向量，含 f1 f2…）
decK           = Dec(idxInSheet, :);     % 该行决策变量

[~, perm] = sort(Obj(:, whichObj), 1, sortDir);
perm = perm(:);

R = struct();
R.file = xlsxPath;
R.sheets = {objSheet, decSheet};
R.nSol = nSol;
R.whichObjective = whichObj;
R.sortDir = sortDir;
R.sortPerm = perm;
R.rank = (1:nSol)';  % 第 k 名
R.Obj = Obj(perm, :);
R.Dec = Dec(perm, :);
R.ObjUnsorted = Obj;
R.DecUnsorted = Dec;
% 说明: R.sortPerm(k) = 全表重排后 第 k 个解 在 原始表中的 行号
R.whichCol   = whichCol;
R.kSmallest  = kSmallest;
R.idxK       = idxInSheet;
R.objK       = objK;
R.decK       = decK;

if readClust && hasSheet('clust_data')
    C = readmatrix(xlsxPath, 'Sheet', 'clust_data', 'NumHeaderLines', 1);
    if size(C, 1) == nSol && size(C, 2) >= 3
        R.clustTag = C(perm, 3);
    elseif size(C, 1) == nSol && size(C, 2) >= 1
        R.clustTag = C(perm, end);
    else
        R.clustTag = [];
        warning('clust_data 尺寸异常，已跳过。');
    end
    R.sheets{end+1} = 'clust_data';
else
    R.clustTag = [];
end

R.details = table();
R.detailsK = table();  % 与「第 whichCol 第 kSmallest 小」同原始行的 details 一行
if readDetails && hasSheet('details_data')
    try
        try
            Td = readtable(xlsxPath, 'Sheet', 'details_data', 'VariableNamingRule', 'preserve');
        catch
            Td = readtable(xlsxPath, 'Sheet', 'details_data');
        end
        if height(Td) == nSol
            R.detailsK = Td(idxInSheet, :);
            R.details = Td(perm, :);
        else
            warning('details_data 行数与 obj_data 不一致，已跳过。');
        end
    catch ME
        warning('读取 details_data 失败: %s', ME.message);
    end
    R.sheets{end+1} = 'details_data';
end

if readHV
    if hasSheet('hv_data')
        R.hv = readmatrix(xlsxPath, 'Sheet', 'hv_data');
    end
    if hasSheet('igd_data')
        R.igd = readmatrix(xlsxPath, 'Sheet', 'igd_data');
    end
end

clear C hasSheet objSheet decSheet

fprintf('已读: %s\n  解数量=%d\n', R.file, R.nSol);
fprintf('【第 %d 列 第 %d 小】原始行号 idxInSheet=%d\n', whichCol, kSmallest, idxInSheet);
fprintf('  该列上该名次取值 = %.6g  （在 sort 升序中第 %d 个）\n', Obj(idxInSheet, whichCol), kSmallest);
fprintf('  objK (整行目标) = '); fprintf('%.6g  ', objK); fprintf('\n');
fprintf('  decK = '); fprintf('%.4g  ', decK); fprintf('\n');
% details_data 中与 obj 同列名的 E/T 分项（需 readDetails=true 且 sheet 行数一致）
detailVarNames = {'E_drone1_J','E_drone2_J','E_drone3_J','T_drone1_s','T_drone2_s','T_drone3_s'};
if height(R.detailsK) == 1
    vn = R.detailsK.Properties.VariableNames;
    for iv = 1:numel(detailVarNames)
        name = detailVarNames{iv};
        if ismember(name, vn)
            v = R.detailsK.(name)(1);
            if isnumeric(v) && isscalar(v)
                fprintf('  %s = %.6g\n', name, v);
            else
                fprintf('  %s = %s\n', name, string(v));
            end
        else
            fprintf('  %s = (无此列；表变量名: %s)\n', name, strjoin(vn, ', '));
        end
    end
elseif readDetails
    fprintf('  (E_drone*_J / T_drone*_s 未输出：details_data 不可用或行数与 obj 不一致)\n');
end
fprintf('\n');

fprintf('全表 R 按 第 %d 列, 方向=%s  重排，sortPerm(前10): ', R.whichObjective, R.sortDir);
disp(R.sortPerm(1:min(10, R.nSol))');
