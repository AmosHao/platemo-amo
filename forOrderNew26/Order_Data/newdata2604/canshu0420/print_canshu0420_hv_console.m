function print_canshu0420_hv_console(rootDir, varargin)
% print_canshu0420_hv_console
% 直接在控制台输出 (pc,er) 参数实验的 HV_end 统计信息（以 baseline 对比）。
%
% 用法：
%   print_canshu0420_hv_console('D:\...\canshu0420')
%   print_canshu0420_hv_console('D:\...\canshu0420','Baseline',[0.5 0.5],'Alpha',0.05,'Adjust','holm')
%
% 输出列：
%   (pc,er) | n | mean | std | delta(Cliff) | p | pHolm | tag
% tag:
%   BASE: baseline
%   +/* : mean > baseline / significant
%   -/* : mean < baseline / significant

    p = inputParser;
    p.addRequired('rootDir', @(s)ischar(s) || isstring(s));
    p.addParameter('Baseline', [0.5 0.5], @(x)isnumeric(x) && numel(x)==2);
    p.addParameter('Alpha', 0.05, @(x)isnumeric(x) && isscalar(x) && x>0 && x<1);
    p.addParameter('Adjust', 'holm', @(s)ischar(s) || isstring(s)); % holm|none
    p.parse(rootDir, varargin{:});

    rootDir = char(p.Results.rootDir);
    basePc = p.Results.Baseline(1);
    baseEr = p.Results.Baseline(2);
    alpha  = p.Results.Alpha;
    adjust = lower(string(p.Results.Adjust));

    % read all summaries
    d = dir(rootDir);
    d = d([d.isdir]);
    d = d(~ismember({d.name},{'.','..'}));
    if isempty(d)
        error('No subfolders under %s', rootDir);
    end

    recs = struct('pc',{},'er',{},'hv',{},'n',{},'mean',{},'std',{});
    for i = 1:numel(d)
        sub = fullfile(rootDir, d(i).name);
        mats = dir(fullfile(sub,'*_runs.mat'));
        if isempty(mats), continue; end
        S = load(fullfile(sub,mats(1).name),'summary');
        if ~isfield(S,'summary') || ~isfield(S.summary,'hv_end_all'), continue; end
        hv = S.summary.hv_end_all(:);
        [pc, er] = parse_pc_er(S.summary, d(i).name, mats(1).name);
        if isnan(pc) || isnan(er), continue; end
        r.pc = pc; r.er = er;
        r.hv = hv;
        r.n = sum(~isnan(hv));
        r.mean = mean(hv,'omitnan');
        r.std  = std(hv,'omitnan');
        recs(end+1) = r; %#ok<AGROW>
    end
    if isempty(recs)
        error('No usable *_runs.mat under %s', rootDir);
    end

    % baseline
    idx0 = find(arrayfun(@(r)abs(r.pc-basePc)<1e-12 && abs(r.er-baseEr)<1e-12, recs), 1);
    if isempty(idx0)
        error('Baseline (pc,er)=(%g,%g) not found.', basePc, baseEr);
    end
    hv0 = recs(idx0).hv;

    % tests vs baseline
    m = numel(recs);
    p_raw = NaN(m,1);
    delta = NaN(m,1);
    for i = 1:m
        if i == idx0, continue; end
        [pval, dlt] = mw_u_test(hv0, recs(i).hv);
        p_raw(i) = pval;
        delta(i) = dlt;
    end

    p_adj = p_raw;
    if adjust == "holm"
        mask = ~isnan(p_raw);
        p_adj(mask) = holm_adjust(p_raw(mask));
    elseif adjust == "none"
        % keep raw
    else
        error('Adjust must be holm or none.');
    end

    sig = false(m,1);
    sig(~isnan(p_adj)) = p_adj(~isnan(p_adj)) < alpha;

    % sort by pc, er
    [~,ord] = sortrows([[recs.pc]', [recs.er]'], [1 2]);
    recs = recs(ord);
    p_raw = p_raw(ord);
    p_adj = p_adj(ord);
    delta = delta(ord);
    sig   = sig(ord);
    idx0  = find(arrayfun(@(r)abs(r.pc-basePc)<1e-12 && abs(r.er-baseEr)<1e-12, recs), 1);

    fprintf('\n=== HV_end stats (baseline pc=%.3g, er=%.3g; alpha=%.3g; adjust=%s) ===\n', ...
        basePc, baseEr, alpha, adjust);
    fprintf('%8s %8s %4s %12s %12s %10s %12s %12s %6s\n', ...
        'pc','er','n','mean','std','delta','p','pAdj','tag');

    for i = 1:m
        if i == idx0
            tag = 'BASE';
            fprintf('%8.3g %8.3g %4d %12.6f %12.6f %10.3f %12s %12s %6s\n', ...
                recs(i).pc, recs(i).er, recs(i).n, recs(i).mean, recs(i).std, 0, '--', '--', tag);
        else
            if recs(i).mean > recs(idx0).mean, sgn = '+'; else, sgn = '-'; end
            if sig(i), tag = [sgn '*']; else, tag = sgn; end
            fprintf('%8.3g %8.3g %4d %12.6f %12.6f %10.3f %12s %12s %6s\n', ...
                recs(i).pc, recs(i).er, recs(i).n, recs(i).mean, recs(i).std, delta(i), ...
                fmt_p(p_raw(i)), fmt_p(p_adj(i)), tag);
        end
    end
    fprintf('tag: + better than baseline, - worse than baseline, * significant.\n\n');
end

% ---------------- helpers (copied from tex generator) ----------------
function s = fmt_p(p)
    if isnan(p), s = '--'; return; end
    if p < 1e-4, s = '<1e-4'; else, s = sprintf('%.4g', p); end
end

function [pc, er] = parse_pc_er(summary, folderName, fileName)
    pc = NaN; er = NaN;
    if isstruct(summary)
        if isfield(summary,'operamo_pc') && ~isempty(summary.operamo_pc), pc = double(summary.operamo_pc); end
        if isfield(summary,'operamo_er') && ~isempty(summary.operamo_er), er = double(summary.operamo_er); end
    end
    if ~isnan(pc) && ~isnan(er), return; end
    token = regexp([folderName ' ' fileName], 'pc([0-9]*\.?[0-9]+)_er([0-9]*\.?[0-9]+)', 'tokens', 'once');
    if ~isempty(token)
        pc = str2double(token{1});
        er = str2double(token{2});
    end
end

function p_adj = holm_adjust(p)
    p = p(:);
    m = numel(p);
    [ps, idx] = sort(p, 'ascend');
    adj = zeros(m,1);
    for i = 1:m
        adj(i) = (m - i + 1) * ps(i);
    end
    for i = 2:m
        adj(i) = max(adj(i), adj(i-1));
    end
    adj = min(adj, 1);
    p_adj = zeros(m,1);
    p_adj(idx) = adj;
end

function [pval, delta] = mw_u_test(x, y)
    x = x(:); y = y(:);
    x = x(~isnan(x)); y = y(~isnan(y));
    nx = numel(x); ny = numel(y);
    if nx < 2 || ny < 2
        pval = NaN; delta = NaN; return;
    end
    delta = cliffs_delta(x, y);

    z = [x; y];
    [zs, ord] = sort(z);
    r = zeros(size(zs));
    i = 1;
    while i <= numel(zs)
        j = i;
        while j < numel(zs) && zs(j+1) == zs(i)
            j = j + 1;
        end
        r(i:j) = (i + j) / 2;
        i = j + 1;
    end
    ranks = zeros(size(z));
    ranks(ord) = r;
    Rx = sum(ranks(1:nx));
    Ux = Rx - nx*(nx+1)/2;
    mu = nx*ny/2;

    tieCounts = diff([0; find(diff(zs)~=0); numel(zs)]);
    T = sum(tieCounts.^3 - tieCounts);
    sigma = sqrt(nx*ny/12 * ((nx+ny+1) - T/((nx+ny)*(nx+ny-1))));
    if sigma == 0
        pval = 1;
        return;
    end
    zstat = (Ux - mu) / sigma;
    pval = 2 * min(normcdf(zstat), 1 - normcdf(zstat));
end

function d = cliffs_delta(x, y)
    nx = numel(x); ny = numel(y);
    gt = 0; lt = 0;
    for i = 1:nx
        gt = gt + sum(x(i) > y);
        lt = lt + sum(x(i) < y);
    end
    d = (gt - lt) / (nx*ny);
end

