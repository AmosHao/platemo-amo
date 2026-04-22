% 直接点“运行”这个脚本即可。
% 通过下面的开关（true/false）选择要画的图：不弹窗、不保存，只显示图窗。

rootDir = 'D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\newdata2604\canshu0420';

% ---------------- 选择要画哪些图（改这里） ----------------
drawHVHeatmap      = true;   % HV 均值/标准差热力图
drawIGDHeatmap     = false;  % IGD 均值/标准差热力图（不需要）
drawRuntimeHeatmap = true;   % 运行时间均值热力图
drawConvergence    = false;   % Top-K 收敛曲线（HV/IGD：均值±std）
drawPareto         = false;   % Top-K 终代 Pareto 前沿散点（代表性 run）

% Top-K：只影响“收敛曲线/帕累托前沿”
topK = 3;

% Top-K 的排序依据：'HV'（越大越好）或 'IGD'（越小越好）
metricPick = 'HV';
% ----------------------------------------------------------

plots = strings(0,1);
if drawHVHeatmap,      plots(end+1) = "HVHeatmap";      end %#ok<AGROW>
if drawIGDHeatmap,     plots(end+1) = "IGDHeatmap";     end %#ok<AGROW>
if drawRuntimeHeatmap, plots(end+1) = "RuntimeHeatmap"; end %#ok<AGROW>
if drawConvergence,    plots(end+1) = "Convergence";    end %#ok<AGROW>
if drawPareto,         plots(end+1) = "Pareto";         end %#ok<AGROW>

plot_canshu0420_param_sweep(rootDir, ...
    'TopK', topK, ...
    'MetricPick', metricPick, ...
    'Plots', plots, ...
    'Save', false, ...   % 不保存
    'Close', false);     % 不关闭（方便查看）

