% 17 条有向 OD：@planner_travel_maxhjd_newobj_fangzhen + @PRMO_fangzhen
% 使用 parfeval/fetchNext 实现“边完成边立刻保存”，避免 worker 内 save 触发 Transparency violation。
%
% 输出：
% - 每条边一份：edge_<from>_to_<to>.mat
% - 额外：meta.mat

function [R, meta] = run_17path_PRMO_fangzhen_parfeval(varargin)
    % 允许“直接点运行此文件”
    if nargin == 0 && nargout == 0
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
        defaultSaveDir = '/home/haichao/Documents/why/why_fromlev/PlatEMO/forPlan26/fangzhenDATA/0430_17_50';

        [R0, meta0] = run_17path_PRMO_fangzhen_parfeval('N',100,'maxFE',60000,'s',20,'z_flight',800, ...
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
    addParameter(p, 'saveDir', '');
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
    dotss = od.dotss; %#ok<NASGU>
    nDepot = size(od.dotss, 1);
    assert(nDepot == 11, 'order_data_fangzhen: dotss 需为 11 行（1~10 点 + 配送中心）');

    % 有向边列表 K×2：每行 [from_id, to_id]。只跑一条边时例如 edges17 = [5, 10];
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
    assert(size(edges17, 2) == 2 && size(edges17, 1) >= 1, 'edges17 须为 K×2（K>=1），每行一条 OD');

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
        end
    end

    if runParallel
        % 强制使用 local（避免 threads worker 不支持 cd 等限制）
        ppool = gcp('nocreate');
        if ~isempty(ppool) && isa(ppool, 'parallel.ThreadPool')
            delete(ppool);
            ppool = [];
        end
        if isempty(gcp('nocreate'))
            nw = min(numWorkers, ne);
            try
                parpool('local', nw);
            catch
                parpool('local');
            end
        end

        poolNow = gcp('nocreate');
        fprintf('[并行模式] parfeval labs=%d，edges=%d\n', poolNow.NumWorkers, ne);

        futures = parallel.FevalFuture.empty(0, 0);
        for e = 1:ne
            a = edges17(e, 1);
            b = edges17(e, 2);
            futures(e) = parfeval( ...
                @edgeTask_PRMO_fangzhen, 1, e, a, b, dotss, z_flight, N, maxFE, s, M);
        end

        completed = 0;
        while completed < ne
            [fidx, out] = fetchNext(futures);
            completed = completed + 1;

            if out.ok
                R(out.e) = out.edge; %#ok<AGROW>
                if ~isempty(saveDir)
                    edge = out.edge; %#ok<NASGU>
                    outName = sprintf('edge_%d_to_%d.mat', out.edge.from_id, out.edge.to_id);
                    outPath = fullfile(saveDir, outName);
                    save(outPath, 'edge', '-v7.3');
                    fprintf('[主线程] 已保存: %s（%d/%d）\n', outPath, completed, ne);
                end
            else
                fprintf('[失败] edge %d/%d (%d->%d) 错误：%s\n', out.e, ne, out.a, out.b, out.err);
            end
        end
    else
        % 串行：边完成就立刻保存
        for e = 1:ne
            a = edges17(e, 1);
            b = edges17(e, 2);
            [out] = edgeTask_PRMO_fangzhen(e, a, b, dotss, z_flight, N, maxFE, s, M); %#ok<ASGLU>
            if out.ok
                R(e) = out.edge;
                if ~isempty(saveDir)
                    edge = out.edge; %#ok<NASGU>
                    outName = sprintf('edge_%d_to_%d.mat', out.edge.from_id, out.edge.to_id);
                    outPath = fullfile(saveDir, outName);
                    save(outPath, 'edge', '-v7.3');
                    fprintf('[主线程] 已保存: %s\n', outPath);
                end
            else
                fprintf('[失败] edge %d/%d (%d->%d) 错误：%s\n', out.e, ne, out.a, out.b, out.err);
            end
        end
    end

    uu = false;
    if isfield(od, 'use_union_routes')
        uu = logical(od.use_union_routes);
    end

    if runParallel
        parallel_backend = 'parfeval+fetchNext(local)';
    else
        parallel_backend = 'serial';
    end

    meta = struct( ...
        'description', 'run_17path_PRMO_fangzhen_parfeval', 'orderFile', orderFile, 'edges17', edges17, ...
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

function out = edgeTask_PRMO_fangzhen(e, a, b, dotss, z_flight, N, maxFE, s, M)
    out.e = e;
    out.a = a;
    out.b = b;
    out.ok = false;
    out.edge = [];
    out.err = '';

    try
        if a == 0, ra = 11; else, ra = a; end
        if b == 0, rb = 11; else, rb = b; end
        p1 = [dotss(ra, 1:2), z_flight];
        p2 = [dotss(rb, 1:2), z_flight];
        dots = [p1, p2];

        [ve, vt, cN, cT, cNe, cTe, Dc, Ob, Dend, Oend, hv, hvE, pdi, pde] = ...
            runplatemo_all_singlepath_PRMO_fangzhen( ...
            @planner_travel_maxhjd_newobj_fangzhen, @PRMO_fangzhen, N, maxFE, s, M, dots);

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

        out.edge = edge;
        out.ok = true;
    catch ME
        out.err = ME.message;
        try
            out.err = getReport(ME, 'extended', 'hyperlinks', 'off');
        catch
        end
    end
end

function od = load_order_data_fangzhen(orderFile)
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

