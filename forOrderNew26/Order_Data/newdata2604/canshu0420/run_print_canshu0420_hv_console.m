% 直接运行本脚本：在控制台输出 HV_end 统计检验结果（vs baseline）

addpath('d:\PlatEMO-master-using\figs');

rootDir = 'D:\PlatEMO-master-using\PlatEMO\forOrderNew26\Order_Data\newdata2604\canshu0420';

% baseline 参数（与你对比试验一致）
baseline = [0.5 0.25];   % [pc, er]

% 显著性水平与多重校正
alpha  = 0.05;
adjust = 'holm';        % 'holm' 或 'none'

print_canshu0420_hv_console(rootDir, 'Baseline', baseline, 'Adjust', adjust, 'Alpha', alpha);

