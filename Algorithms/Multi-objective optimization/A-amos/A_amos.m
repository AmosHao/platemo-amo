classdef A_amos < ALGORITHM
% <multi> <real/integer>
% A-amos: Adaptive Multi-objective Optimization for Order Assignment
% 基于PRMO改造的订单分配多目标优化算法
% 
% 参数说明:
% Kmax --- 20 --- 目标空间聚类最大数量（随种群增大而增加）
% BETA --- 0.8 --- 邻域选择平衡参数
% pc --- 0.9 --- 交叉概率（保持高交叉率以充分利用大种群）
% pm --- 0.2 --- 变异概率（适度增加以提高探索能力）
% maxgen --- 300 --- 最大迭代次数

    methods
        function main(Algorithm, Problem)
           % 屏蔽复数绘图警告
           warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
           
           %% 参数设置
           % 调整参数以适应更大的种群：
           % - Kmax: 聚类数随种群增大而增加
           % - pc/pm: 保持较高的交叉概率，增加变异概率以提高多样性
           % - maxgen: 增加迭代次数以充分利用大种群
           [Kmax, BETA, pc, pm, maxgen] = Algorithm.ParameterSet(20, 0.8, 0.8, 0.3, 300);
           
           % 修复概率（降低以增加多样性，避免过度修复导致局部最优）
           repairProb = 0.6;  % 60%概率进行顺序约束修复
           completenessProb = 0.7;  % 70%概率进行完整性修复
           
           popSize = Problem.N; 
           objDim = Problem.M;
           varDim = Problem.D;
           n = Problem.n;  % 客户点数量（商家数量）
           m = Problem.m;  % 无人机数量
           
           %% 初始化种群
           Population = Problem.Initialization();
           pop = Population.decs;
           objVals = Population.objs;
           
           % 初始化聚类信息（在目标空间）
           clustTag = (1:popSize)';
           clustName = (1:popSize)';
           centroid = objVals;  % 聚类中心在目标空间
           
           %% 局部最优避免机制：停滞检测与自适应参数
           hvHistory = [];          % 记录历史超体积（HV）值，越大越好
           stagnationCounter = 0;   % 停滞计数器
           stagnationThreshold = 40;% 停滞阈值（连续40代HV无显著增长，温和化 30–50）
           adaptivePm = pm;         % 自适应变异概率
           adaptivePc = pc;         % 自适应交叉概率
           restartInterval = 50;    % 重启间隔（每50代检查一次）
           diversityHistory = [];   % 多样性历史
           
           %% 精英档案：永久保存历史最优端点解，防止好的端点被截断淘汰
           eliteF1Dec = [];   % f1目标最小的解（决策向量）
           eliteF1Obj = [];   % f1目标最小的解（目标值）
           eliteF2Dec = [];   % f2目标最小的解（决策向量）
           eliteF2Obj = [];   % f2目标最小的解（目标值）
           
               %% 主优化循环
           gen = 0;
           while Algorithm.NotTerminated(Population)
               gen = gen + 1;
               
              %% 停滞检测与自适应参数调整（基于超体积HV）— 已注释，禁用自适应停滞
              % if ~isempty(objVals) && size(objVals, 1) > 0
              %     validObj = objVals(objVals(:,1) < 1e9, :);
              %     if ~isempty(validObj)
              %         currentHV = compute_hv_2d(validObj);
              %         hvHistory = [hvHistory, currentHV];
              %         if size(validObj, 1) > 1
              %             diversity = mean(std(validObj, 0, 1)) / (mean(mean(validObj)) + 1e-10);
              %             diversityHistory = [diversityHistory, diversity];
              %         end
              %         if length(hvHistory) > stagnationThreshold
              %             recentHV = hvHistory(end-stagnationThreshold+1:end);
              %             improvement = (recentHV(end) - recentHV(1)) / (abs(recentHV(1)) + 1e-10);
              %             if improvement < 0.025
              %                 stagnationCounter = stagnationCounter + 1;
              %                 adaptivePm = min(0.5, pm * (1 + stagnationCounter * 0.08));
              %                 adaptivePc = min(0.95, pc * (1 + stagnationCounter * 0.02));
              %                 if mod(gen, 20) == 0
              %                     fprintf('[Gen %3d] 停滞(计数=%d): HV=%.4f, 增长=%.4f%%, Pm=%.3f, Pc=%.3f\n', ...
              %                         gen, stagnationCounter, currentHV, improvement*100, adaptivePm, adaptivePc);
              %                 end
              %             else
              %                 stagnationCounter = max(0, stagnationCounter - 1);
              %                 adaptivePm = pm;
              %                 adaptivePc = pc;
              %             end
              %         end
              %         if length(diversityHistory) > 10
              %             recentDiversity = mean(diversityHistory(end-9:end));
              %             if recentDiversity < 0.01
              %                 adaptivePm = min(0.5, adaptivePm * 1.1);
              %                 adaptivePc = min(0.95, adaptivePc * 1.05);
              %                 if mod(gen, 20) == 0
              %                     fprintf('[Gen %3d] 多样性过低: diversity=%.6f, Pm=%.3f, Pc=%.3f\n', ...
              %                         gen, recentDiversity, adaptivePm, adaptivePc);
              %                 end
              %             end
              %         end
              %     end
              % end
               
               %% 随机重启机制（恢复有效性：一次性替换10%，打破局部最优）
               % 注意：重启时的新解会通过Initialization评估，但这是必要的
               % 主循环中只评估auxPop（子代），保证每代只调用一次Evaluation（在EnvironmentalSelection中）
               % if mod(gen, restartInterval) == 0 && stagnationCounter > 5
               %     % 替换10%的最差解为随机新解（批量，保持速度）
               %     if ~isempty(objVals) && size(objVals, 1) > 0
               %         numReplace = max(1, floor(popSize * 0.1));
               %         [~, worstIdx] = sort(objVals(:,1), 'descend');
               %         replaceIdx = worstIdx(1:numReplace);
               %         newPop = Problem.Initialization(numReplace);
               %         pop(replaceIdx, :) = newPop.decs;
               %         objVals(replaceIdx, :) = newPop.objs;
               %         fprintf('代 %d: 随机重启，替换 %d 个最差解（已评估）\n', gen, numReplace);
               %         stagnationCounter = 0;  % 重置停滞计数器
               %     end
               % end
               
               %% 目标空间聚类管理
               auxPop = [];
               globalClust = [];  % 获取全局子种群
               nClust = size(clustName, 1);  % 聚类数量
               neigSet = cell(nClust, 1);
               
               % 确保pop和objVals不为空
               if isempty(pop) || isempty(objVals)
                   continue;  % 如果种群为空，跳过本次迭代
               end
               
               % 为每个聚类构建邻域集合（基于目标空间距离）
               for i = 1:nClust
                   if i > length(clustName)
                       continue;  % 跳过无效索引
                   end
                   mark = clustTag == clustName(i);  % 找到当前属于第i类的个体编号
                   if sum(mark) > 0
                       neigPop = pop(mark, :);
                       neigObj = objVals(mark, :);
                       
                       % 从每个聚类取 f1 最小和 f2 最小的解加入全局池（若同一解则只加一次）
                       [~, posF1] = min(neigObj(:, 1));
                       [~, posF2] = min(neigObj(:, 2));
                       globalClust = [globalClust; neigPop(posF1, :)];
                       if posF2 ~= posF1
                           globalClust = [globalClust; neigPop(posF2, :)];
                       end
                       
                       if sum(mark) > 2
                           neigSet{i, 1} = {neigPop, neigObj};  % 存储邻域解和目标值
                       end
                   end
               end
               % 方法一：全局父代从全种群选，保留聚类仅用于邻域内选，避免全局池过小、多人同父代
               globalClust = pop;
               globalSize = size(globalClust, 1);
               
             %% 个体进化
             % [算子调试计数器] 每代重置
             dbg_cross_ptc  = 0;  % 位置追踪交叉(PTC)调用次数
             dbg_cross_ppc  = 0;  % 对优先级交叉(PPC)调用次数
             dbg_cross_opc  = 0;  % 顺序保持交叉(OPC)调用次数
             dbg_cross_def  = 0;  % 交叉触发克隆防御次数
             dbg_mut_t1     = 0;  % 变异策略1调用次数（商家前移）
             dbg_mut_t2     = 0;  % 变异策略2调用次数（客户后移）
             dbg_mut_t3     = 0;  % 变异策略3调用次数（分隔符平移）
             dbg_mut_t4     = 0;  % 变异策略4调用次数（跨段对迁移）
             dbg_mut_def    = 0;  % 变异触发克隆防御次数
             for i = 1:popSize
                   if i > size(pop, 1) || i > size(objVals, 1) || i > length(clustTag)
                       continue;  % 跳过无效索引
                   end
                   currentSol = pop(i, :);
                   currentObj = objVals(i, :);
                   currentTag = clustTag(i);
                   
                   % 确保currentSol维度正确
                   if length(currentSol) ~= varDim
                       % if length(currentSol) < varDim
                       %     currentSol = [currentSol, zeros(1, varDim - length(currentSol))];
                       % else
                       %     currentSol = currentSol(1:varDim);
                       % end
                       fprintf('Warning: 维度不正确');
                       
                   end
                   
                   % 获取当前个体的邻域
                   neighborhood = [];
                   neighborhoodObj = [];
                   clustIdx = find(clustName == currentTag, 1);
                   if ~isempty(clustIdx) && clustIdx <= length(neigSet) && ~isempty(neigSet{clustIdx, 1})
                       temp = neigSet{clustIdx, 1};
                       neighborhood = temp{1};
                       neighborhoodObj = temp{2};
                       % 从邻域中删除当前个体
                       if ~isempty(neighborhood)
                           [~, loc] = ismember(currentSol, neighborhood, 'rows');
                           if loc > 0 && loc <= size(neighborhood, 1)
                               neighborhood(loc, :) = [];
                               neighborhoodObj(loc, :) = [];
                           end
                       end
                   end
                   neigSize = size(neighborhood, 1);
                   
                  %% 父代逻辑：交叉=当前+不同类随机一个，变异=当前个体（自身）；按概率二选一，固定 0.5
                  % [已注释] 按迭代调节概率：genRatio = gen/maxgen; crossoverProb = max(0, pc*(1-genRatio)); ...
                  % genRatio       = gen / maxgen;
                  % crossoverProb  = max(0, pc * (1 - genRatio));
                  % mutationProb   = max(0.4, 1 - crossoverProb);
                  crossoverProb  = 0.5;   % 固定 0.5
                  useCrossover   = rand < crossoverProb;
                  % useMutation    = rand < mutationProb;
                  explorationRate = min(0.70, 0.25 + stagnationCounter * 0.02);
                  
                 % 交叉：当前个体 + 不同类的另一个；变异：当前类中随机一个
                 if useCrossover && popSize >= 2
                     parent = currentSol;
                     otherClustIdx = find(clustTag ~= currentTag);
                     otherClustIdx = setdiff(otherClustIdx, i);
                     if isempty(otherClustIdx)
                         otherClustIdx = setdiff(1:popSize, i);
                     end
                     idx2 = otherClustIdx(randsample(length(otherClustIdx), 1));
                     parent2 = pop(idx2, :);
                     if isequal(parent2, currentSol) && length(otherClustIdx) > 1
                         idx2 = otherClustIdx(randsample(length(otherClustIdx), 1));
                         parent2 = pop(idx2, :);
                     end
                     [trialSol, crossDbg] = SegmentBasedCrossover(parent, parent2, n, m, varDim);
                     if crossDbg.strategy == 1, dbg_cross_ptc = dbg_cross_ptc + 1;
                     elseif crossDbg.strategy == 2, dbg_cross_ppc = dbg_cross_ppc + 1;
                     else, dbg_cross_opc = dbg_cross_opc + 1; end
                     if crossDbg.cloneDefended, dbg_cross_def = dbg_cross_def + 1; end
                     if length(trialSol) ~= varDim
                         trialSol = [trialSol, zeros(1, max(0, varDim - length(trialSol)))];
                         trialSol = trialSol(1:varDim);
                     end
                 else
                     % 变异：父代固定为当前个体（自身）
                     parent = currentSol;
                     [trialSol, mutDbg] = OrderPreservingMutation(parent, m, varDim, n, explorationRate);
                     if mutDbg.strategy == 1, dbg_mut_t1 = dbg_mut_t1 + 1;
                     elseif mutDbg.strategy == 2, dbg_mut_t2 = dbg_mut_t2 + 1;
                     elseif mutDbg.strategy == 3, dbg_mut_t3 = dbg_mut_t3 + 1;
                     else, dbg_mut_t4 = dbg_mut_t4 + 1; end
                     if mutDbg.cloneDefended, dbg_mut_def = dbg_mut_def + 1; end
                 end
                   
                  % 增强探索性操作：避免局部最优
                  % 1. 当停滞时增加随机扰动（段交换）- 仅当两段均为完整对（偶数长度）时交换，避免产生 "0 0 18"
                  % if stagnationCounter > 2 && rand < 0.4
                  %     idx0 = find(trialSol == 0);
                  %     if length(idx0) >= 2
                  %         segIdx = randperm(length(idx0)+1, 2);
                  %         idx0_exp = [0, idx0, length(trialSol)+1];
                  %         seg1_start = idx0_exp(segIdx(1)) + 1;
                  %         seg1_end = idx0_exp(segIdx(1)+1) - 1;
                  %         seg2_start = idx0_exp(segIdx(2)) + 1;
                  %         seg2_end = idx0_exp(segIdx(2)+1) - 1;
                  %         if seg1_end >= seg1_start && seg2_end >= seg2_start
                  %             seg1 = trialSol(seg1_start:seg1_end);
                  %             seg2 = trialSol(seg2_start:seg2_end);
                  %             % 仅当两段长度相等且均为偶数（完整商家-客户对）时交换
                  %             L1 = length(seg1); L2 = length(seg2);
                  %             if L1 == L2 && mod(L1, 2) == 0
                  %                 trialSol(seg1_start:seg1_end) = seg2;
                  %                 trialSol(seg2_start:seg2_end) = seg1;
                  %             end
                  %         end
                  %     end
                  % end
                % 2. 定期（每10代）对所有个体增加一次额外变异，打破局部最优（主算法中关闭该额外变异）
                % if mod(gen, 10) == 0 && rand < 0.2
                %     [trialSol, mutDbg2] = OrderPreservingMutation(trialSol, m, varDim, n, explorationRate);
                %      if mutDbg2.strategy == 1, dbg_mut_t1 = dbg_mut_t1 + 1;
                %      elseif mutDbg2.strategy == 2, dbg_mut_t2 = dbg_mut_t2 + 1;
                %      elseif mutDbg2.strategy == 3, dbg_mut_t3 = dbg_mut_t3 + 1;
                %      else, dbg_mut_t4 = dbg_mut_t4 + 1; end
                %      if mutDbg2.cloneDefended, dbg_mut_def = dbg_mut_def + 1; end
                % end
                   
                   % 统一段结构修复：确保每段均为完整商家-客户对（可两 0 相邻，但左右段必须完整）
                   trialSol = OrderSegmentRepair(trialSol, n, m, varDim);
                   
                   % 确保是行向量
                   if size(trialSol, 1) ~= 1
                       trialSol = reshape(trialSol, 1, []);
                   end
                   
                   % 添加到auxPop
                   if isempty(auxPop)
                       auxPop = trialSol;
                   else
                       auxPop = [auxPop; trialSol];
                   end
               end
               
             %% [算子调试] 每10代输出一次算子调用统计
             if mod(gen, 10) == 0
                 crossTotal = dbg_cross_ptc + dbg_cross_ppc + dbg_cross_opc;
                 mutTotal   = dbg_mut_t1 + dbg_mut_t2 + dbg_mut_t3 + dbg_mut_t4;
                 fprintf('[Gen %3d] 交叉(%d次): PTC=%d PPC=%d OPC=%d | 克隆防御=%d\n', ...
                     gen, crossTotal, dbg_cross_ptc, dbg_cross_ppc, dbg_cross_opc, dbg_cross_def);
                 fprintf('[Gen %3d] 变异(%d次):  T1=%d  T2=%d  T3=%d  T4=%d | 克隆防御=%d  explorationRate=%.2f\n', ...
                     gen, mutTotal, dbg_mut_t1, dbg_mut_t2, dbg_mut_t3, dbg_mut_t4, dbg_mut_def, explorationRate);
             end
              
              %% 环境选择前：缓存精英对应的 Population 对象（在种群变化前取，零额外FE）
              % 精英解已在历史迭代中评估，目标值存于 eliteF1Obj/eliteF2Obj，不需重新评估
              eliteF1PopObj = [];
              eliteF2PopObj = [];
              if ~isempty(eliteF1Dec)
                  f1InPop = find(all(pop == eliteF1Dec, 2), 1);
                  if ~isempty(f1InPop) && f1InPop <= length(Population)
                      eliteF1PopObj = Population(f1InPop);
                  end
              end
              if ~isempty(eliteF2Dec) && ~isequal(eliteF2Dec, eliteF1Dec)
                  f2InPop = find(all(pop == eliteF2Dec, 2), 1);
                  if ~isempty(f2InPop) && f2InPop <= length(Population)
                      eliteF2PopObj = Population(f2InPop);
                  end
              end
              
              %% 环境选择：只传入 auxPop（N个子代），每代固定消耗 2N 次 FE
              selectedSize = Problem.N;
              [pop, objVals, clustTag, clustName, centroid, Population] = ...
                  EnvironmentalSelection(auxPop, pop, objVals, clustTag, clustName, ...
                  centroid, selectedSize, objDim, varDim, Kmax, Problem);
              
              popSize = size(pop, 1);
              
              %% 更新精英档案（从选择后种群中提取当代最优端点，单调更新）
              if ~isempty(objVals)
                  validMask = objVals(:,1) < 1e9;
                  if any(validMask)
                      validObj = objVals(validMask, :);
                      validDec = pop(validMask, :);
                      [~, f1MinIdx] = min(validObj(:,1));
                      if isempty(eliteF1Obj) || validObj(f1MinIdx,1) < eliteF1Obj(1)
                          eliteF1Dec = validDec(f1MinIdx, :);
                          eliteF1Obj = validObj(f1MinIdx, :);
                          % 同步更新缓存的 Population 对象
                          f1InNew = find(all(pop == eliteF1Dec, 2), 1);
                          if ~isempty(f1InNew) && f1InNew <= length(Population)
                              eliteF1PopObj = Population(f1InNew);
                          end
                      end
                      [~, f2MinIdx] = min(validObj(:,2));
                      if isempty(eliteF2Obj) || validObj(f2MinIdx,2) < eliteF2Obj(2)
                          eliteF2Dec = validDec(f2MinIdx, :);
                          eliteF2Obj = validObj(f2MinIdx, :);
                          f2InNew = find(all(pop == eliteF2Dec, 2), 1);
                          if ~isempty(f2InNew) && f2InNew <= length(Population)
                              eliteF2PopObj = Population(f2InNew);
                          end
                      end
                  end
              end
              
              %% 强制精英保留：若精英在环境选择后被淘汰，写回种群替换最差解，零额外FE
              if ~isempty(eliteF1Dec) && ~isempty(objVals) && popSize > 0
                  if ~any(all(pop == eliteF1Dec, 2))
                      [~, worstIdx] = max(objVals(:,1));
                      pop(worstIdx, :) = eliteF1Dec;
                      objVals(worstIdx, :) = eliteF1Obj;
                      clustTag(worstIdx) = inf;
                      if ~isempty(eliteF1PopObj) && worstIdx <= length(Population)
                          Population(worstIdx) = eliteF1PopObj;
                      end
                  end
              end
              if ~isempty(eliteF2Dec) && ~isempty(objVals) && popSize > 0 && ~isequal(eliteF2Dec, eliteF1Dec)
                  if ~any(all(pop == eliteF2Dec, 2))
                      [~, worstIdx] = max(objVals(:,2));
                      pop(worstIdx, :) = eliteF2Dec;
                      objVals(worstIdx, :) = eliteF2Obj;
                      clustTag(worstIdx) = inf;
                      if ~isempty(eliteF2PopObj) && worstIdx <= length(Population)
                          Population(worstIdx) = eliteF2PopObj;
                      end
                  end
              end
               
               % 更新Population对象
               if popSize > 0 && length(Population) >= popSize
                   for i = 1:popSize
                       if i <= length(Population)
                           Population(i).add_clustTag = clustTag(i, :);
                           Population(i).add_clustName = clustName(:, :);
                       end
                   end
               end
           end
        end
    end
end

function hv = compute_hv_2d(objVals)
% 计算二维超体积（用于停滞检测）
% 使用固定参考范围归一化，与 HV.m 保持一致，但先过滤惩罚无效解
% 输入：objVals - N×2 目标值矩阵（调用前已过滤 f1<1e9）
% 输出：hv - 归一化后的超体积值（[0,1] 范围）
    fmin = [0, 0];
    fmax = [1e7, 4e5];  % 与 HV.m 中 fmax 保持一致

    % 归一化到 [0, 1/1.1] 空间
    N = size(objVals, 1);
    normObj = (objVals - repmat(fmin, N, 1)) ./ repmat((fmax - fmin) * 1.1, N, 1);
    % 只剔除负值（不应出现的情况），对超出上界的解裁剪到边界而非直接删除
    % 这样避免因少量解轻微越界而导致 HV 突然跌落
    normObj(any(normObj < 0, 2), :) = [];
    if isempty(normObj), hv = 0; return; end
    normObj = min(normObj, 1.0);  % 裁剪到参考点边界，不删除

    % 提取二维非支配前沿（O(n log n)）
    normObj = sortrows(normObj, [1, 2]);  % 先按 f1 升序
    n = size(normObj, 1);
    keep = true(n, 1);
    min_f2 = inf;
    for i = 1:n
        if normObj(i, 2) < min_f2
            min_f2 = normObj(i, 2);
        else
            keep(i) = false;  % 被前面的点支配
        end
    end
    nd = normObj(keep, :);  % 非支配前沿，已按 f1 升序

    % 扫描线算法计算 2D HV（参考点 [1,1]）
    % 每个非支配点贡献一个水平条带：
    %   宽度 = RefPoint(1) - nd(i,1)
    %   高度 = prev_f2 - nd(i,2)
    hv = 0;
    prev_f2 = 1.0;  % 参考点 f2 = 1.0
    for i = 1:size(nd, 1)
        if nd(i, 2) < prev_f2
            hv = hv + (1.0 - nd(i, 1)) * (prev_f2 - nd(i, 2));
            prev_f2 = nd(i, 2);
        end
    end
end

function completeSol = EnsureSolutionCompleteness(sol, n, m, varDim)
% 确保解包含所有必需的点（1到2*n）
% **关键修复**：确保0分隔符数量正确（m-1个），并在调整维度时优先保留0分隔符

    completeSol = sol;
    
    % 必需的点：1到2*n（商家1-10，客户点11-20）
    allRequiredPoints = 1:(2*n);
    
    % 找到所有非0的点
    nonZeroVals = completeSol(completeSol ~= 0);
    uniqueVals = unique(nonZeroVals);
    
    % 检查是否缺少必需的点
    missingVals = setdiff(allRequiredPoints, uniqueVals);
    
    if ~isempty(missingVals)
        % 如果缺少点，需要补充
        % 策略：替换重复的点，或插入到0分隔符附近
        
        % 找出重复的点
        if ~isempty(nonZeroVals)
            [counts, ~, idx] = unique(nonZeroVals);
            freq = histcounts(idx, 1:max(idx)+1);
            duplicateVals = nonZeroVals(freq > 1);
        else
            duplicateVals = [];
        end
        
        % 先替换重复的点
        dupIdx = 1;
        for i = 1:length(completeSol)
            if completeSol(i) ~= 0 && ismember(completeSol(i), duplicateVals)
                % 检查是否已经出现过
                prevOccurrences = find(completeSol(1:i-1) == completeSol(i));
                if ~isempty(prevOccurrences) && dupIdx <= length(missingVals)
                    completeSol(i) = missingVals(dupIdx);
                    dupIdx = dupIdx + 1;
                end
            end
        end
        
        % 如果还有缺失的点，插入到0分隔符附近
        remainingMissing = setdiff(allRequiredPoints, unique(completeSol(completeSol ~= 0)));
        if ~isempty(remainingMissing)
            zeroIdx = find(completeSol == 0);
            if ~isempty(zeroIdx) && zeroIdx(1) > 1
                % 在第一个0之前插入
                completeSol = [completeSol(1:zeroIdx(1)-1), remainingMissing, completeSol(zeroIdx(1):end)];
            else
                % 在末尾插入
                completeSol = [completeSol, remainingMissing];
            end
        end
    end
    
    % **关键修复**：先确保0分隔符数量正确（m-1个），再调整维度
    % 这样可以避免在调整维度时丢失0分隔符
    zeroCount = sum(completeSol == 0);
    if zeroCount < m - 1
        % 如果0分隔符不足，需要插入
        % 策略：在非0元素之间均匀插入0分隔符
        nonZeroIdx = find(completeSol ~= 0);
        if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
            % 计算需要插入的0的数量
            numZerosNeeded = (m - 1) - zeroCount;
            % 在非0元素之间均匀插入0
            insertPositions = round(linspace(1, length(nonZeroIdx)-1, numZerosNeeded+1));
            insertPositions = insertPositions(2:end);  % 去掉第一个位置
            % 从后往前插入，避免索引变化
            for i = length(insertPositions):-1:1
                pos = nonZeroIdx(insertPositions(i));
                completeSol = [completeSol(1:pos), 0, completeSol(pos+1:end)];
            end
        else
            % 如果没有非0元素或只有一个，直接在末尾添加0
            completeSol = [completeSol, zeros(1, (m-1) - zeroCount)];
        end
    elseif zeroCount > m - 1
        % 如果0分隔符过多，只保留前m-1个
        zeroIdx = find(completeSol == 0);
        if length(zeroIdx) > m - 1
            % 删除多余的0，但要保留前m-1个
            toRemove = zeroIdx(m:end);
            completeSol(toRemove) = [];
        end
    end
    
    % 确保维度正确（在确保0分隔符之后）
    % **关键**：调整维度时，要确保0分隔符不被截断
    currentZeroCount = sum(completeSol == 0);
    if length(completeSol) ~= varDim
        if length(completeSol) < varDim
            % 如果长度不足，补充0（但要确保0分隔符数量不超过m-1）
            remainingZeros = (m-1) - currentZeroCount;
            if remainingZeros > 0
                % 先补充0分隔符
                completeSol = [completeSol, zeros(1, remainingZeros)];
            end
            % 再补充到目标维度（用非0值填充，但这里我们用0，因为varDim应该已经包含了0分隔符）
            % 实际上，varDim应该等于 2*n + (m-1)，所以如果长度不足，应该是缺少非0元素
            % 但为了安全，我们只补充到varDim，不再添加额外的0
            if length(completeSol) < varDim
                completeSol = [completeSol, zeros(1, varDim - length(completeSol))];
            end
        else
            % 如果长度超过varDim，需要截断
            % **关键**：截断时要确保保留m-1个0分隔符
            zeroIdx = find(completeSol == 0);
            if length(zeroIdx) >= m - 1
                % 保留前m-1个0分隔符，截断其他部分
                lastZeroIdx = zeroIdx(m-1);
                % 计算需要保留的非0元素数量
                nonZeroNeeded = varDim - (m - 1);
                nonZeroCount = sum(completeSol(1:lastZeroIdx) ~= 0);
                if nonZeroCount <= nonZeroNeeded
                    % 保留到最后一个0分隔符之后，再补充一些非0元素
                    completeSol = completeSol(1:min(varDim, length(completeSol)));
                else
                    % 需要截断一些非0元素
                    completeSol = completeSol(1:varDim);
                end
            else
                % 如果0分隔符不足，先补充0，再截断
                completeSol = completeSol(1:varDim);
            end
        end
    end
    
    % 最终验证：确保0分隔符数量正确
    finalZeroCount = sum(completeSol == 0);
    if finalZeroCount ~= m - 1
        % 如果还是不对，强制修正
        zeroIdx = find(completeSol == 0);
        if finalZeroCount < m - 1
            % 补充0分隔符
            numNeeded = (m - 1) - finalZeroCount;
            % 在非0元素之间插入
            nonZeroIdx = find(completeSol ~= 0);
            if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
                insertPos = round(linspace(1, length(nonZeroIdx)-1, numNeeded+1));
                insertPos = insertPos(2:end);
                for i = length(insertPos):-1:1
                    pos = nonZeroIdx(insertPos(i));
                    completeSol = [completeSol(1:pos), 0, completeSol(pos+1:end)];
                end
            else
                completeSol = [completeSol, zeros(1, numNeeded)];
            end
        elseif finalZeroCount > m - 1
            % 删除多余的0
            zeroIdx = find(completeSol == 0);
            completeSol(zeroIdx(m:end)) = [];
        end
        % 再次调整维度
        if length(completeSol) ~= varDim
            if length(completeSol) < varDim
                completeSol = [completeSol, zeros(1, varDim - length(completeSol))];
            else
                completeSol = completeSol(1:varDim);
            end
        end
    end
end
