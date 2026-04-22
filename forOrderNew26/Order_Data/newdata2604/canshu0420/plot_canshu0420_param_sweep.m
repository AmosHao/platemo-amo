% 运行
% plot_canshu0420_param_sweep('D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\newdata2604\canshu0420')

function plot_canshu0420_param_sweep(rootDir, varargin)
% plot_canshu0420_param_sweep
% 一键汇总 NSGAIII_G 的 (pc,er) 参数实验结果，并生成论文常用图：
% - HV_end / IGD_end 的 mean/std 热力图（pc×er）
% - 运行时间 mean 热力图
% - 参数组的收敛曲线（均值±标准差）
% - Top-K 参数组的一次代表性终代 Pareto 前沿散点（选 HV 最好的 run）
%
% 用法：
%   plot_canshu0420_param_sweep('D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\newdata2604\canshu0420')
%
% 可选参数（Name-Value）：
%   'TopK'        : 默认 3
%   'MetricPick'  : 'HV' 或 'IGD'（决定 Top-K 的排序依据，默认 'HV'）
%   'OutDir'      : 输出图片目录（默认 rootDir 下 fig_param_sweep）
%   'Plots'       : 要画的图（cell 或 string 数组），默认 "all"
%                  可选：
%                    - "HVHeatmap"      (HV 均值/标准差热力图)
%                    - "IGDHeatmap"     (IGD 均值/标准差热力图)
%                    - "RuntimeHeatmap" (运行时间均值热力图)
%                    - "Convergence"    (全部参数组收敛曲线：HV 均值±std)
%                    - "Pareto"         (Top-K 终代前沿散点)
%   'Save'        : 是否保存为 png（默认 true）
%   'Close'       : 保存后是否关闭图窗（默认 true）。Save=false 时建议 Close=false 以便查看
%
% 说明：
% - 脚本优先使用 summary.operamo_pc / summary.operamo_er；若缺失则从文件夹/文件名解析。
% - 假设每个参数组目录下只有 1 个 *_runs.mat。

    p = inputParser;
    p.addRequired('rootDir', @(s)ischar(s) || isstring(s));
    p.addParameter('TopK', 3, @(x)isnumeric(x) && isscalar(x) && x>=1);
    p.addParameter('MetricPick', 'HV', @(s)ischar(s) || isstring(s));
    p.addParameter('OutDir', '', @(s)ischar(s) || isstring(s));
    p.addParameter('Plots', "all", @(x)ischar(x) || isstring(x) || iscell(x));
    p.addParameter('Save', true, @(x)islogical(x) && isscalar(x));
    p.addParameter('Close', true, @(x)islogical(x) && isscalar(x));
    p.parse(rootDir, varargin{:});
    topK = p.Results.TopK;
    metricPick = upper(string(p.Results.MetricPick));
    if metricPick ~= "HV" && metricPick ~= "IGD"
        error('MetricPick must be ''HV'' or ''IGD''.');
    end

    plots = normalize_plots(p.Results.Plots);
    doSave = p.Results.Save;
    doClose = p.Results.Close;

    rootDir = char(rootDir);
    if isempty(p.Results.OutDir)
        outDir = fullfile(rootDir, 'fig_param_sweep');
    else
        outDir = char(p.Results.OutDir);
    end
    if doSave
        if ~isfolder(outDir), mkdir(outDir); end
    end

    % ---- enumerate all *_runs.mat in subfolders ----
    d = dir(rootDir);
    d = d([d.isdir]);
    d = d(~ismember({d.name},{'.','..'}));
    if isempty(d)
        error('No subfolders found under: %s', rootDir);
    end

    records = struct('pc',{},'er',{},'matPath',{},'hv_end_all',{},'igd_end_all',{}, ...
        'hv_curves',{},'igd_curves',{},'run_time_sec_all',{},'obj_all_runs',{},'runs',{});

    for i = 1:numel(d)
        subDir = fullfile(rootDir, d(i).name);
        mats = dir(fullfile(subDir, '*_runs.mat'));
        if isempty(mats)
            % 兼容：有时 mat 放在 rootDir 直接一层
            continue;
        end
        if numel(mats) > 1
            warning('Multiple *_runs.mat found in %s, using first: %s', subDir, mats(1).name);
        end
        matPath = fullfile(subDir, mats(1).name);
        S = load(matPath, 'summary', 'runs');
        if ~isfield(S,'summary')
            warning('Missing summary in %s, skipped.', matPath);
            continue;
        end
        summary = S.summary;

        % pc/er
        [pc, er] = get_pc_er(summary, d(i).name, mats(1).name);

        rec = struct();
        rec.pc = pc;
        rec.er = er;
        rec.matPath = matPath;
        rec.hv_end_all = getfield_def(summary, 'hv_end_all', []);
        rec.igd_end_all = getfield_def(summary, 'igd_end_all', []);
        rec.hv_curves = getfield_def(summary, 'hv_curves', []);
        rec.igd_curves = getfield_def(summary, 'igd_curves', []);
        rec.run_time_sec_all = getfield_def(summary, 'run_time_sec_all', []);
        rec.obj_all_runs = getfield_def(summary, 'obj_all_runs', []);
        if isfield(S,'runs'), rec.runs = S.runs; else, rec.runs = []; end

        records(end+1) = rec; %#ok<AGROW>
    end

    if isempty(records)
        error('No *_runs.mat found under: %s', rootDir);
    end

    pcs = unique([records.pc]);
    ers = unique([records.er]);
    pcs = sort(pcs);
    ers = sort(ers);

    % ---- build grid stats ----
    HV_mean = NaN(numel(pcs), numel(ers));
    HV_std  = NaN(numel(pcs), numel(ers));
    IGD_mean = NaN(numel(pcs), numel(ers));
    IGD_std  = NaN(numel(pcs), numel(ers));
    T_mean  = NaN(numel(pcs), numel(ers));

    for k = 1:numel(records)
        ip = find(pcs == records(k).pc, 1);
        ie = find(ers == records(k).er, 1);
        if ~isempty(records(k).hv_end_all)
            HV_mean(ip,ie) = mean(records(k).hv_end_all(:),'omitnan');
            HV_std(ip,ie)  = std(records(k).hv_end_all(:),'omitnan');
        end
        if ~isempty(records(k).igd_end_all)
            IGD_mean(ip,ie) = mean(records(k).igd_end_all(:),'omitnan');
            IGD_std(ip,ie)  = std(records(k).igd_end_all(:),'omitnan');
        end
        if ~isempty(records(k).run_time_sec_all)
            T_mean(ip,ie) = mean(records(k).run_time_sec_all(:),'omitnan');
        end
    end

    % ---- heatmaps ----
    if want_plot(plots, "HVHeatmap")
        save_heatmap(HV_mean, pcs, ers, 'HV均值', fullfile(outDir,'HV_mean.png'), doSave, doClose);
        save_heatmap(HV_std,  pcs, ers, 'HV标准差',  fullfile(outDir,'HV_std.png'),  doSave, doClose);
    end
    if want_plot(plots, "IGDHeatmap")
        save_heatmap(IGD_mean, pcs, ers, 'IGD\_end mean', fullfile(outDir,'IGD_mean.png'), doSave, doClose);
        save_heatmap(IGD_std,  pcs, ers, 'IGD\_end std',  fullfile(outDir,'IGD_std.png'),  doSave, doClose);
    end
    if want_plot(plots, "RuntimeHeatmap")
        save_heatmap(T_mean,  pcs, ers, '运行时间 (s)', fullfile(outDir,'runtime_mean.png'), doSave, doClose);
    end

    % ---- pick Top-K parameter groups (only if needed) ----
    needTopK = want_plot(plots, "Pareto");
    if needTopK
        score = NaN(1, numel(records));
        for k = 1:numel(records)
            if metricPick == "HV"
                if ~isempty(records(k).hv_end_all)
                    score(k) = mean(records(k).hv_end_all(:),'omitnan');
                end
            else
                if ~isempty(records(k).igd_end_all)
                    score(k) = mean(records(k).igd_end_all(:),'omitnan');
                end
            end
        end
        valid = find(~isnan(score));
        if isempty(valid)
            warning('No valid metric data found for Top-K plots.');
            if doSave
                fprintf('[OK] Figures saved to: %s\n', outDir);
            else
                fprintf('[OK] Figures shown (not saved).\n');
            end
            return;
        end
        if metricPick == "HV"
            [~,ord] = sort(score(valid), 'descend');
        else
            [~,ord] = sort(score(valid), 'ascend');
        end
        chosen = valid(ord(1:min(topK,numel(ord))));
    else
        chosen = [];
    end

    % ---- convergence curves (ALL parameter groups; HV only) ----
    if want_plot(plots, "Convergence")
        f = figure('Color','w','Position',[100 100 760 520]);
        hold on; grid on;
        title('HV convergence for all (pc,er) (mean \pm std)');
        xlabel('Generation');
        ylabel('HV');
        for idx = 1:numel(records)
            [mu, sd] = curve_mean_std(records(idx).hv_curves);
            if isempty(mu), continue; end
            plot_with_band(mu, sd, sprintf('pc=%.3g, er=%.3g', records(idx).pc, records(idx).er));
        end
        legend('Location','eastoutside');
        if doSave
            exportgraphics(f, fullfile(outDir, 'convergence_all.png'), 'Resolution', 300);
        end
        if doClose
            close(f);
        end
    end

    % ---- Pareto fronts (representative run) ----
    if want_plot(plots, "Pareto")
        f = figure('Color','w','Position',[100 100 520 520]);
        hold on; grid on; box on;
        title(sprintf('Final Pareto fronts (Top-%d by %s)', numel(chosen), metricPick));
        xlabel('f1');
        ylabel('f2');
        colors = lines(numel(chosen));
        for j = 1:numel(chosen)
            idx = chosen(j);
            [Obj, label] = representative_front(records(idx), metricPick);
            if isempty(Obj), continue; end
            scatter(Obj(:,1), Obj(:,2), 12, 'MarkerEdgeColor', colors(j,:), 'DisplayName', label);
        end
        legend('Location','best');
        if doSave
            exportgraphics(f, fullfile(outDir, 'pareto_topk.png'), 'Resolution', 300);
        end
        if doClose
            close(f);
        end
    end

    if doSave
        fprintf('[OK] Figures saved to: %s\n', outDir);
    else
        fprintf('[OK] Figures shown (not saved).\n');
    end
end

% ---------------- helpers ----------------
function plots = normalize_plots(x)
    % x can be char/string/cell -> normalize to string array
    plots = string(x);
    plots = upper(strtrim(plots));
    if any(plots == "ALL")
        plots = ["HVHEATMAP","IGDHEATMAP","RUNTIMEHEATMAP","CONVERGENCE","PARETO"];
    end
end

function tf = want_plot(plots, key)
    tf = any(plots == upper(string(key)));
end

function v = getfield_def(s, name, def)
    if isstruct(s) && isfield(s, name)
        v = s.(name);
    else
        v = def;
    end
end

function [pc, er] = get_pc_er(summary, folderName, fileName)
    pc = NaN; er = NaN;
    if isstruct(summary)
        if isfield(summary,'operamo_pc') && ~isempty(summary.operamo_pc), pc = summary.operamo_pc; end
        if isfield(summary,'operamo_er') && ~isempty(summary.operamo_er), er = summary.operamo_er; end
    end
    if ~isnan(pc) && ~isnan(er), return; end

    % parse from strings like "..._pc0.25_er0.75..."
    token = regexp([folderName ' ' fileName], 'pc([0-9]*\.?[0-9]+)_er([0-9]*\.?[0-9]+)', 'tokens', 'once');
    if isempty(token)
        token = regexp([folderName ' ' fileName], 'pc([0-9]*\.?[0-9]+)\s*[_-]\s*er([0-9]*\.?[0-9]+)', 'tokens', 'once');
    end
    if ~isempty(token)
        pc = str2double(token{1});
        er = str2double(token{2});
    end
end

function save_heatmap(A, pcs, ers, ttl, outPath, doSave, doClose)
    f = figure('Color','w','Position',[100 100 560 480]);
    % 用均匀的下标 1..n 画 n×n 色块，避免 imagesc(数值轴) 与 grid 对单元格进行“对半切割”
    ny = size(A,1);
    nx = size(A,2);
    hi = imagesc(1:nx, 1:ny, A);
    set(gca,'YDir','normal');
    colormap(parula);
    colorbar;
    xlabel('$\rho$','Interpreter','latex');
    ylabel('$P_c$','Interpreter','latex');
    title(ttl,'Interpreter','tex');
    xticks(1:nx);
    yticks(1:ny);
    if exist('string','class')
        xticklabels(string(ers(:)'));
        yticklabels(string(pcs(:)'));
    else
        set(gca,'XTickLabel',arrayfun(@num2str, ers, 'UniformOutput', false));
        set(gca,'YTickLabel',arrayfun(@num2str, pcs, 'UniformOutput', false));
    end
    axis([0.5, nx+0.5, 0.5, ny+0.5]);
    set(gca,'Box','on','LineWidth',1.5);
    % 网格画在相邻格子的分界处 (k+0.5)，不要对 grid on（格线落在刻度/格心把色块切开）
    hold on;
    xl = [0.5, nx+0.5];
    yl = [0.5, ny+0.5];
    for xv = 1.5:1:(nx-0.5)
        plot([xv xv], yl, ':', 'Color', [0.92 0.92 0.92], 'LineWidth', 0.9, 'Clipping', 'on', 'HandleVisibility', 'off');
    end
    for yv = 1.5:1:(ny-0.5)
        plot(xl, [yv yv], ':', 'Color', [0.92 0.92 0.92], 'LineWidth', 0.9, 'Clipping', 'on', 'HandleVisibility', 'off');
    end
    hold off;
    if ~isempty(hi) && isgraphics(hi)
        uistack(hi, 'bottom');
    end
    if nargin < 6, doSave = true; end
    if nargin < 7, doClose = true; end
    if doSave
        exportgraphics(f, outPath, 'Resolution', 300);
    end
    if doClose
        close(f);
    end
end

function [mu, sd] = curve_mean_std(curves)
    mu = []; sd = [];
    if isempty(curves), return; end
    % curves is (nGen x nRuns), may contain NaN padding
    mu = mean(curves, 2, 'omitnan');
    sd = std(curves, 0, 2, 'omitnan');
end

function plot_with_band(mu, sd, label)
    x = (1:numel(mu))';
    y1 = mu - sd;
    y2 = mu + sd;
    h = plot(x, mu, 'LineWidth', 1.6, 'DisplayName', label);
    c = h.Color;
    fill([x; flipud(x)], [y1; flipud(y2)], c, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility','off');
end

function [Obj, label] = representative_front(rec, metricPick)
% 选一个代表性 run 的前沿：HV 取最大；IGD 取最小
    Obj = [];
    label = sprintf('pc=%.3g, er=%.3g', rec.pc, rec.er);
    if nargin < 2 || isempty(metricPick), metricPick = "HV"; end
    metricPick = upper(string(metricPick));

    if isempty(rec.runs)
        if ~isempty(rec.obj_all_runs) && size(rec.obj_all_runs,2) >= 3
            Obj = rec.obj_all_runs(:,2:3);
            label = sprintf('%s (all runs)', label);
        end
        return;
    end
    if metricPick == "IGD"
        vals = arrayfun(@(s)s.igd_end, rec.runs);
        [~,best] = min(vals);
        label = sprintf('%s (bestIGD run %d)', label, best);
    else
        vals = arrayfun(@(s)s.hv_end, rec.runs);
        [~,best] = max(vals);
        label = sprintf('%s (bestHV run %d)', label, best);
    end
    Obj = rec.runs(best).Obj_end;
end

