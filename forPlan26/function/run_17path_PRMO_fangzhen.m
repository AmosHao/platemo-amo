% 17 条有向 OD：@planner_travel_maxhjd_newobj_fangzhen + @PRMO_fangzhen
% 支持并行（推荐 spmd+laboratory；默认开）：参考 ParTest_all_TRAVEL_PRMO.m -> SpmdRun_all_TRAVEL_PRMO_class.m
% 每条边在独立 lab 上跑，结果写入各自 edge_*.mat
% 输出结构与 runplatemo_all_singlepath_PRMO_class / ParTest 类似
%
% 节点编号（与 order_data_fangzhen 一致）：0 = 配送中心（dotss 第 11 行），1~10 = dotss 第 1~10 行
% 航迹端点高度使用 z_flight（与问题类 z_min 一致，默认 800 m）

function [R, meta] = run_17path_PRMO_fangzhen(varargin)
    % 允许“直接点运行此文件”：无参无输出时，自动按默认参数执行并把结果写回工作区
    if nargin == 0 && nargout == 0
        % 不写死路径：从本文件位置推导工程根目录
        thisFile = mfilename('fullpath');
        ocDir = fileparts(thisFile);
        forOrderDir = fileparts(ocDir);
        plDir = fileparts(forOrderDir);
        repoRoot = fileparts(plDir);
        if isempty(repoRoot)
            repoRoot = pwd;
        end
        plRoot = fullfile(repoRoot, 'PlatEMO');
        addpath(plRoot, genpath(plRoot));
        setenv('PLATEMO_FANGZHEN_ROOT', repoRoot);

        % 你要改“保存到哪里”，就改这一行
        % defaultSaveDir = fullfile(repoRoot, 'OUT', 'PRMO17');
        defaultSaveDir = '/home/haichao/Documents/why/why_fromlev/PlatEMO/forPlan26/fangzhenDATA/0429_17';

        [R0, meta0] = run_17path_PRMO_fangzhen('N',100,'maxFE',600,'s',20,'z_flight',800, ...
            'saveDir', defaultSaveDir);
        assignin('base','R',R0);
        assignin('base','meta',meta0);
        fprintf('已在工作区生成变量：R, meta\n');
        return;
    end

    p = inputParser;
    addParameter(p, 'N', 100);
    addParameter(p, 'maxFE', 60000);
    addParameter(p, 's', 20);
    addParameter(p, 'M', 3);
    addParameter(p, 'z_flight', 800);
    % 若 saveDir 非空，则输出 17 个文件：edge_<from>_to_<to>.mat，以及 meta.mat
    addParameter(p, 'saveDir', '');
    % 并行：每条边一任务；默认开启。使用 local 进程池 + spmd（与仓库内 SpmdRun_* 脚本一致）。
    addParameter(p, 'useParallel', true);
    addParameter(p, 'numWorkers', 17);
    parse(p, varargin{:});
    N = p.Results.N;
    maxFE = p.Results.maxFE;
    s = p.Results.s;
    M = p.Results.M;
    z_flight = p.Results.z_flight;
    saveDir = p.Results.saveDir;
    useParallel = p.Results.useParallel;
    numWorkers = p.Results.numWorkers;

    thisFile = mfilename('fullpath');
    ocDir = fileparts(thisFile);
    forOrderDir = fileparts(ocDir);
    plDir = fileparts(forOrderDir);
    repoRoot = fileparts(plDir);
    if isempty(repoRoot)
        repoRoot = pwd;
    end
    setenv('PLATEMO_FANGZHEN_ROOT', repoRoot);

    plRoot = fullfile(repoRoot, 'PlatEMO');
    if ~isfile(fullfile(plRoot, 'platemo.m'))
        error('未找到 %s，请从工程根运行或检查路径', fullfile(plRoot, 'platemo.m'));
    end
    addpath(plRoot, genpath(plRoot));
    addpath(fullfile(plRoot, 'forOrderNew26', 'Order_Map'));
    addpath(fullfile(plRoot, 'forPlan26', 'function'));

    orderFile = fullfile(plRoot, 'forOrderNew26', 'Order_Map', 'order_data_fangzhen.m');
    od = load_order_data_fangzhen(orderFile);
    dotss = od.dotss;
    nDepot = size(dotss, 1);
    assert(nDepot == 11, 'order_data_fangzhen: dotss 需为 11 行（1~10 点 + 配送中心）');

    %% 17 条有向边（两条路径合并去重后的 17 个有序对；D=配送中心=0）
    % 路 A 独有：D->1, 6->2, 7->4, 4->9
    % 路 B 独有：D->4, 4->2, 7->1, 6->9
    % 两条都用到：D->3, 3->8, 8->D, D->5, 5->10, 10->D, 2->7, 1->6, 9->D
    % 顺序可按需要调整：这里把 4->9 放到第一条先跑
    edges17 = [
        4, 9;
        0, 1;
        6, 2;
        7, 4;
        0, 4;
        4, 2;
        7, 1;
        6, 9;
        0, 3;
        3, 8;
        8, 0;
        0, 5;
        5, 10;
        10, 0;
        2, 7;
        1, 6;
        9, 0
    ];
    assert(size(edges17, 1) == 17, '当前 edges17 行数 = %d，需改为 17 对', size(edges17, 1));

    R = repmat(struct( ...
        'from_id', 0, 'to_id', 0, 'dots', [], ...
        'validtrait_end', [], 'validtrait', [], 'clustName', [], 'clustTag', [], ...
        'clustName_end', [], 'clustTag_end', [], ...
        'Dec', [], 'Obj', [], 'Dec_end', [], 'Obj_end', [], ...
        'hv', [], 'hv_end', NaN, 'pd', [], 'pd_end', NaN), size(edges17, 1), 1);

    if ~isempty(saveDir)
        if ~isfolder(saveDir)
            mkdir(saveDir);
        end
    end

    ne = size(edges17, 1);
    runParallel = useParallel;
    if runParallel
        if exist('parpool', 'file') ~= 2 || ~license('test', 'Distrib_Computing_Toolbox')
            warning('未安装或未授权 Parallel Computing Toolbox，改为串行。');
            runParallel = false;
        else
            try
                ppool = gcp('nocreate');
                nw = min(numWorkers, ne);
                % 这里必须优先使用 local 进程池：线程池 worker 不支持 cd，而 PlatEMO 调用链里常会 cd
                %（如 platemo / 初始化脚本）。若当前已存在线程池，则切换为 local。
                if ~isempty(ppool) && isa(ppool, 'parallel.ThreadPool')
                    delete(ppool);
                    ppool = [];
                end
                if isempty(ppool)
                    try
                        parpool('local', nw);
                    catch
                        % 本机允许的最大 lab 数可能小于 nw（如仅 4/8 核许可）
                        parpool('local');
                    end
                end

                % 若是 Processes 池（进程），给 worker 补齐环境与路径；Threads 池与主线程共享，无需处理
                poolNow = gcp('nocreate');
                isThreadsPool = false;
                if ~isempty(poolNow)
                    % 兼容不同 MATLAB 版本：ThreadPool 可能没有 Type 属性
                    if isa(poolNow, 'parallel.ThreadPool')
                        isThreadsPool = true;
                    elseif isprop(poolNow, 'Type')
                        isThreadsPool = strcmpi(poolNow.Type, 'Threads');
                    end
                end
                if ~isempty(poolNow) && ~isThreadsPool
                    rr = repoRoot;
                    plR = plRoot;
                    cmd1 = sprintf('setenv(''PLATEMO_FANGZHEN_ROOT'',''%s'');', rr);
                    cmd2 = sprintf('addpath(''%s'', genpath(''%s''));', plR, plR);
                    cmd3 = sprintf('addpath(''%s'');', fullfile(plR, 'forOrderNew26', 'Order_Map'));
                    cmd4 = sprintf('addpath(''%s'');', fullfile(plR, 'forPlan26', 'function'));
                    if exist('pctRunOnAll', 'file') == 2
                        pctRunOnAll(cmd1);
                        pctRunOnAll(cmd2);
                        pctRunOnAll(cmd3);
                        pctRunOnAll(cmd4);
                    elseif exist('parfevalOnAll', 'file') == 2
                        % 避免嵌套 spmd：用 parfevalOnAll 在所有 worker 执行初始化脚本串
                        parfevalOnAll(@eval, 0, [cmd1 cmd2 cmd3 cmd4]); %#ok<PFENV>
                    end
                end
            catch ME
                warning('失败: %s，改为串行。', ME.message);
                runParallel = false;
            end
        end
    end

    if runParallel
        poolNow = gcp('nocreate');
        actualLabs = poolNow.NumWorkers;
        fprintf('[并行模式] spmd labs=%d（请求 workers=%d），edges=%d\n', actualLabs, min(numWorkers, ne), ne);

        % 不显式传入 labs 数量：自动使用当前 parallel pool 的 worker 数（与 SpmdRun_all_TRAVEL_PRMO_class 一致）
        spmd
            localEdges = [];
            for e = labindex:numlabs:ne
                tEdge = tic;
                rng(e + 1e4 + labindex, 'twister');
                a = edges17(e, 1);
                b = edges17(e, 2);
                stage = 'prepare_inputs';
                try
                    if a == 0, ra = 11; else, ra = a; end
                    if b == 0, rb = 11; else, rb = b; end
                    p1 = [dotss(ra, 1:2), z_flight];
                    p2 = [dotss(rb, 1:2), z_flight];
                    dots = [p1, p2];

                    fprintf('[lab %d/%d] 开始 %d -> %d\n', labindex, numlabs, a, b);

                    stage = 'solve_edge';
                    [ve, vt, cN, cT, cNe, cTe, Dc, Ob, Dend, Oend, hv, hvE, pdi, pde] = ...
                        runplatemo_all_singlepath_PRMO_fangzhen( ...
                        @planner_travel_maxhjd_newobj_fangzhen, @PRMO_fangzhen, N, maxFE, s, M, dots);

                    stage = 'pack_edge_struct';
                    edge = struct();
                    edge.from_id = a;
                    edge.to_id = b;
                    edge.dots = dots;
                    edge.validtrait_end = ve;
                    edge.validtrait = vt;
                    edge.clustName = cN;
                    edge.clustTag = cT;
                    edge.clustName_end = cNe;
                    edge.clustTag_end = cTe;
                    edge.Dec = Dc;
                    edge.Obj = Ob;
                    edge.Dec_end = Dend;
                    edge.Obj_end = Oend;
                    edge.hv = hv;
                    edge.hv_end = hvE;
                    edge.pd = pdi;
                    edge.pd_end = pde;

                    % 注意：worker 内 save() 会触发 Transparency violation（已验证）。
                    % 这里仅计算并回传 edge，统一在主线程落盘。
                    fprintf('[lab %d/%d] 完成 %d -> %d | %.1fs\n', labindex, numlabs, a, b, toc(tEdge));

                    localEdges = [localEdges; struct('idx', e, 'edge', edge)]; %#ok<AGROW>
                catch MEe
                    errStage = stage;
                    errReport = MEe.message;
                    try
                        errReport = getReport(MEe, 'extended', 'hyperlinks', 'off');
                    catch
                    end
                    fprintf('[lab %d/%d] 失败 %d -> %d | %.1fs | 阶段=%s | 错误: %s\n', labindex, numlabs, a, b, toc(tEdge), errStage, errReport);
                    bad = struct();
                    bad.idx = e;
                    bad.edge = [];
                    bad.err = MEe;
                    localEdges = [localEdges; bad]; %#ok<AGROW>
                end
            end
        end

        % 组装回主线程：Composite localEdges -> 串起来（SpmdRun_* 同款思路）
        merged = [];
        for lab = 1:actualLabs
            le = localEdges{lab};
            if ~isempty(le)
                merged = [merged; le]; %#ok<AGROW>
            end
        end

        for k = 1:numel(merged)
            item = merged(k);
            if isempty(item.edge)
                continue;
            end
            ee = item.idx;
            R(ee).from_id = item.edge.from_id;
            R(ee).to_id = item.edge.to_id;
            R(ee).dots = item.edge.dots;
            R(ee).validtrait_end = item.edge.validtrait_end;
            R(ee).validtrait = item.edge.validtrait;
            R(ee).clustName = item.edge.clustName;
            R(ee).clustTag = item.edge.clustTag;
            R(ee).clustName_end = item.edge.clustName_end;
            R(ee).clustTag_end = item.edge.clustTag_end;
            R(ee).Dec = item.edge.Dec;
            R(ee).Obj = item.edge.Obj;
            R(ee).Dec_end = item.edge.Dec_end;
            R(ee).Obj_end = item.edge.Obj_end;
            R(ee).hv = item.edge.hv;
            R(ee).hv_end = item.edge.hv_end;
            R(ee).pd = item.edge.pd;
            R(ee).pd_end = item.edge.pd_end;
        end

        % 统一在主线程落盘（避免 worker 的 save() 透明性限制）
        if ~isempty(saveDir)
            for e = 1:ne
                edge = R(e); %#ok<NASGU>
                outName = sprintf('edge_%d_to_%d.mat', R(e).from_id, R(e).to_id);
                outPath = fullfile(saveDir, outName);
                save(outPath, 'edge', '-v7.3');
                fprintf('[主线程] 已保存: %s\n', outPath);
            end
        end
        clear merged

    else
        for e = 1:ne
            a = edges17(e, 1);
            b = edges17(e, 2);
            if a == 0, ra = 11; else, ra = a; end
            if b == 0, rb = 11; else, rb = b; end
            p1 = [dotss(ra, 1:2), z_flight];
            p2 = [dotss(rb, 1:2), z_flight];
            dots = [p1, p2];
            fprintf('==== 第 %d/%d: %d -> %d (dots 1x6) ====\n', e, ne, a, b);

            [ve, vt, cN, cT, cNe, cTe, Dc, Ob, Dend, Oend, hv, hvE, pdi, pde] = ...
                runplatemo_all_singlepath_PRMO_fangzhen( ...
                @planner_travel_maxhjd_newobj_fangzhen, @PRMO_fangzhen, N, maxFE, s, M, dots);

            R(e).from_id = a;
            R(e).to_id = b;
            R(e).dots = dots;
            R(e).validtrait_end = ve;
            R(e).validtrait = vt;
            R(e).clustName = cN;
            R(e).clustTag = cT;
            R(e).clustName_end = cNe;
            R(e).clustTag_end = cTe;
            R(e).Dec = Dc;
            R(e).Obj = Ob;
            R(e).Dec_end = Dend;
            R(e).Obj_end = Oend;
            R(e).hv = hv;
            R(e).hv_end = hvE;
            R(e).pd = pdi;
            R(e).pd_end = pde;

            if ~isempty(saveDir)
                edge = R(e); %#ok<NASGU>
                outName = sprintf('edge_%d_to_%d.mat', a, b);
                outPath = fullfile(saveDir, outName);
                save(outPath, 'edge', '-v7.3');
                fprintf('已保存: %s\n', outPath);
            end
        end
    end

    uu = false;
    if isfield(od, 'use_union_routes')
        uu = logical(od.use_union_routes);
    end
    if runParallel
        parallel_backend = 'spmd+gcp(local)';
    else
        parallel_backend = 'serial';
    end
    meta = struct( ...
        'description', 'run_17path_PRMO_fangzhen', 'orderFile', orderFile, 'edges17', edges17, ...
        'N', N, 'maxFE', maxFE, 's', s, 'M', M, 'z_flight', z_flight, ...
        'use_parallel', runParallel, 'num_workers_requested', numWorkers, ...
        'parallel_backend', parallel_backend, ...
        'use_union_routes', uu);

    if ~isempty(saveDir)
        metaPath = fullfile(saveDir, 'meta.mat');
        save(metaPath, 'meta', 'repoRoot', '-v7.3');
        fprintf('已保存: %s\n', metaPath);
    end

end

function od = load_order_data_fangzhen(orderFile)
% 在函数工作区执行脚本，并把需要的变量导出为结构体。
% 这样 parfor 能稳定识别 dotss 为普通变量（避免 eval 产生的静态/不可分析工作区问题）。
    % 注意：工程内可能存在自定义 run.m 覆盖 MATLAB 内置 run，因此这里不用 run()。
    % 直接读文件并在本函数工作区 eval，避免路径覆盖问题。
    txt = fileread(orderFile);
    eval(txt); %#ok<EVLDIR>
    if ~exist('dotss', 'var')
        error('order_data_fangzhen 未在脚本中生成变量 dotss：%s', orderFile);
    end
    od = struct();
    od.dotss = dotss;
    if exist('use_union_routes', 'var')
        od.use_union_routes = use_union_routes;
    end
end
