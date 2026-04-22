function tex = gen_hv_significance_table_tex(rootDir, pc0, er0, varargin)
% gen_hv_significance_table_tex
% 读取参数实验目录下所有 *_runs.mat，以 (pc0,er0) 为基线，
% 对 HV_end 做两两显著性检验（Mann–Whitney U, two-sided）并输出 booktabs 三线表 LaTeX。
%
% 用法：
%   tex = gen_hv_significance_table_tex(rootDir, 0.5, 0.5)
%   tex = gen_hv_significance_table_tex(rootDir, 0.5, 0.5, 'Alpha',0.05, 'Adjust','holm')
%
% 输出：
%   tex : char，LaTeX 表格源码（booktabs 风格）
%
% 说明：
% - 每个参数组目录下应有 1 个 *_runs.mat（包含 summary.hv_end_all）
% - p 值默认做 Holm 校正（可选 'none'）

    p = inputParser;
    p.addRequired('rootDir', @(s)ischar(s) || isstring(s));
    p.addRequired('pc0', @(x)isnumeric(x) && isscalar(x));
    p.addRequired('er0', @(x)isnumeric(x) && isscalar(x));
    p.addParameter('Alpha', 0.05, @(x)isnumeric(x) && isscalar(x) && x>0 && x<1);
    p.addParameter('Adjust', 'holm', @(s)ischar(s) || isstring(s)); % holm|none
    p.addParameter('Caption', 'Significance test on $HV_{end}$ (vs. baseline).', @(s)ischar(s) || isstring(s));
    p.addParameter('Label', 'tab:hv_sig', @(s)ischar(s) || isstring(s));
    p.parse(rootDir, pc0, er0, varargin{:});

    rootDir = char(p.Results.rootDir);
    alpha = p.Results.Alpha;
    adjust = lower(string(p.Results.Adjust));
    cap = char(p.Results.Caption);
    lab = char(p.Results.Label);
    cap = sanitize_caption(cap);

    % enumerate subfolders
    d = dir(rootDir);
    d = d([d.isdir]);
    d = d(~ismember({d.name},{'.','..'}));
    if isempty(d), error('No subfolders under %s', rootDir); end

    recs = struct('pc',{},'er',{},'hv',{},'mean',{},'std',{},'n',{},'path',{});
    for i = 1:numel(d)
        sub = fullfile(rootDir, d(i).name);
        mats = dir(fullfile(sub,'*_runs.mat'));
        if isempty(mats), continue; end
        matPath = fullfile(sub, mats(1).name);
        S = load(matPath,'summary');
        if ~isfield(S,'summary') || ~isfield(S.summary,'hv_end_all'), continue; end
        hv = S.summary.hv_end_all(:);
        [pc,er] = parse_pc_er(S.summary, d(i).name, mats(1).name);
        if isnan(pc) || isnan(er), continue; end
        r = struct();
        r.pc = pc; r.er = er; r.hv = hv;
        r.mean = mean(hv,'omitnan');
        r.std = std(hv,'omitnan');
        r.n = sum(~isnan(hv));
        r.path = matPath;
        recs(end+1) = r; %#ok<AGROW>
    end
    if isempty(recs), error('No usable *_runs.mat found under %s', rootDir); end

    % baseline
    idx0 = find(arrayfun(@(r)abs(r.pc-pc0)<1e-12 && abs(r.er-er0)<1e-12, recs), 1);
    if isempty(idx0)
        error('Baseline (pc,er)=(%g,%g) not found under %s', pc0, er0, rootDir);
    end
    hv0 = recs(idx0).hv;

    % tests
    m = numel(recs);
    p_raw = NaN(m,1);
    delta = NaN(m,1);
    dirc  = strings(m,1);
    for i = 1:m
        if i == idx0
            p_raw(i) = NaN;
            delta(i) = 0;
            dirc(i)  = "=";
            continue;
        end
        [pval, dlt] = mw_u_test(hv0, recs(i).hv);
        p_raw(i) = pval;
        delta(i) = dlt;
        if dlt > 0, dirc(i) = "+"; elseif dlt < 0, dirc(i) = "-"; else, dirc(i) = "="; end
    end

    % adjust p
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

    % sort rows by pc then er
    [~,ord] = sortrows([[recs.pc]', [recs.er]'], [1 2]);
    recs = recs(ord);
    p_raw = p_raw(ord);
    p_adj = p_adj(ord);
    delta = delta(ord);
    dirc  = dirc(ord);
    sig   = sig(ord);
    idx0_new = find(arrayfun(@(r)abs(r.pc-pc0)<1e-12 && abs(r.er-er0)<1e-12, recs), 1);

    % build tex
    lines = strings(0,1);
    lines(end+1) = "\\begin{table}[t]";
    lines(end+1) = "\\centering";
    lines(end+1) = sprintf("\\caption{%s}", cap);
    lines(end+1) = sprintf("\\label{%s}", lab);
    lines(end+1) = "\\small";
    lines(end+1) = "\\begin{tabular}{ccccccl}";
    lines(end+1) = "\\toprule";
    if adjust == "holm"
        lines(end+1) = "$(pc,er)$ & $n$ & $\\mu\\pm\\sigma$ & $p$ & $p_{Holm}$ & $\\delta$ & Result \\\\";
    else
        lines(end+1) = "$(pc,er)$ & $n$ & $\\mu\\pm\\sigma$ & $p$ & $\\delta$ & Result \\\\";
    end
    lines(end+1) = "\\midrule";

    for i = 1:m
        isBase = (i == idx0_new);
        muSig = sprintf("%.6f$\\pm$%.6f", recs(i).mean, recs(i).std);
        if isBase
            if adjust == "holm"
                row = sprintf("(%.3g, %.3g) & %d & %s & -- & -- & 0.000 & Baseline \\\\", ...
                    recs(i).pc, recs(i).er, recs(i).n, muSig);
            else
                row = sprintf("(%.3g, %.3g) & %d & %s & -- & 0.000 & Baseline \\\\", ...
                    recs(i).pc, recs(i).er, recs(i).n, muSig);
            end
        else
            if isnan(p_raw(i))
                prow = "--";
            else
                prow = format_p(p_raw(i));
            end
            if adjust == "holm"
                padj = format_p(p_adj(i));
                tag = sprintf("%s%s", dirc(i), ternary(sig(i),"*",""));
                row = sprintf("(%.3g, %.3g) & %d & %s & %s & %s & %.3f & %s \\\\", ...
                    recs(i).pc, recs(i).er, recs(i).n, muSig, prow, padj, delta(i), tag);
            else
                tag = sprintf("%s%s", dirc(i), ternary(p_raw(i)<alpha,"*",""));
                row = sprintf("(%.3g, %.3g) & %d & %s & %s & %.3f & %s \\\\", ...
                    recs(i).pc, recs(i).er, recs(i).n, muSig, prow, delta(i), tag);
            end
        end
        lines(end+1) = row;
    end

    lines(end+1) = "\\bottomrule";
    lines(end+1) = "\\end{tabular}";
    if adjust == "holm"
        lines(end+1) = sprintf("\\\\[-1mm]{\\footnotesize $p$ from Mann--Whitney U (two-sided). $p_{Holm}$: Holm correction. $\\delta$: Cliff's delta (vs. baseline). Result: $+$ better / $-$ worse / $=$ tie; $^*$ means significant at $\\alpha=%.2f$.}", alpha);
    else
        lines(end+1) = sprintf("\\\\[-1mm]{\\footnotesize $p$ from Mann--Whitney U (two-sided). $\\delta$: Cliff's delta (vs. baseline). Result: $+$ better / $-$ worse / $=$ tie; $^*$ means significant at $\\alpha=%.2f$.}", alpha);
    end
    lines(end+1) = "\\end{table}";

    tex = strjoin(lines, newline);
end

% ------------------ helpers ------------------
function cap = sanitize_caption(cap)
% Make caption LaTeX-safe:
% - Convert common math tokens like HV_{end} -> $HV_{end}$
% - Escape remaining underscores (outside $...$)
    cap = char(string(cap));

    % Prefer simple literal replacements (robust across MATLAB regex variants)
    cap = strrep(cap, 'HV_{end}', '$HV_{end}$');
    cap = strrep(cap, 'IGD_{end}', '$IGD_{end}$');

    % Protect math segments $...$ before escaping underscores
    mathTokens = {};
    k = 0;
    while true
        [s,e] = regexp(cap, '\$[^$]*\$', 'start', 'end', 'once');
        if isempty(s), break; end
        k = k + 1;
        key = sprintf('<<MATH_%d>>', k);
        mathTokens{k,1} = key; %#ok<AGROW>
        mathTokens{k,2} = cap(s:e); %#ok<AGROW>
        cap = [cap(1:s-1), key, cap(e+1:end)];
    end

    % Escape leftover underscores not already escaped
    cap = regexprep(cap, '(?<!\\)_', '\_');

    % Restore math segments
    for i = 1:size(mathTokens,1)
        cap = strrep(cap, mathTokens{i,1}, mathTokens{i,2});
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function pstr = format_p(p)
    if p < 1e-4
        pstr = "$<10^{-4}$";
    else
        pstr = sprintf("%.4f", p);
    end
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
    % enforce monotonicity
    for i = 2:m
        adj(i) = max(adj(i), adj(i-1));
    end
    adj = min(adj, 1);
    p_adj = zeros(m,1);
    p_adj(idx) = adj;
end

function [pval, delta] = mw_u_test(x, y)
% Mann–Whitney U test (two-sided) + Cliff's delta
% 不依赖 Statistics Toolbox。
    x = x(:); y = y(:);
    x = x(~isnan(x)); y = y(~isnan(y));
    nx = numel(x); ny = numel(y);
    if nx < 2 || ny < 2
        pval = NaN; delta = NaN; return;
    end

    % Cliff's delta
    delta = cliffs_delta(x, y);

    % ranks with ties
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
    % tie correction for variance
    tieCounts = diff([0; find(diff(zs)~=0); numel(zs)]);
    T = sum(tieCounts.^3 - tieCounts);
    sigma = sqrt(nx*ny/12 * ((nx+ny+1) - T/((nx+ny)*(nx+ny-1))));
    if sigma == 0
        pval = 1;
        return;
    end
    zstat = (Ux - mu) / sigma;
    % two-sided using normal approximation
    pval = 2 * min(normcdf(zstat), 1 - normcdf(zstat));
end

function d = cliffs_delta(x, y)
    nx = numel(x); ny = numel(y);
    % O(n^2) is fine for ~31 runs
    gt = 0; lt = 0;
    for i = 1:nx
        gt = gt + sum(x(i) > y);
        lt = lt + sum(x(i) < y);
    end
    d = (gt - lt) / (nx*ny);
end

