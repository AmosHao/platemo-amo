classdef GOCEAfororder_1207< ALGORITHM
% <multi> <real/integer>
% GOCEA
% Kmax ---   11 --- The K
% BETA    --- 0.8 --- The B
% F ---   0.5 --- The F
% CR    --- 1 --- The CR
    methods
        function main(Algorithm, Problem)
           %% parameter setting
           [Kmax,BETA,F,CR] = Algorithm.ParameterSet(11,0.8,0.5,1);
           popSize = Problem.N; 
           objDim=Problem.M;
           varDim=Problem.D;
           % pm=1.0/varDim;
           % bounds=[Problem.lower;Problem.upper];
           %% Generate random population
           Population = Problem.Initialization();
           
           % 多样性跟踪变量
           diversity_history = [];
           stagnation_counter = 0;
           last_best_hv = 0;
           
           %% Optimization
           gen = 1;
           while Algorithm.NotTerminated(Population)
               gen = gen + 1;
               
               % 每一代从 Population 同步 pop / objvs，保证三者一致
               pop   = Population.decs;
               objvs = Population.objs;
               popSize = size(pop,1);
               
               % 多样性检查和自适应调整
               current_diversity = calculate_diversity(pop, objvs);
               diversity_history = [diversity_history, current_diversity];
               
               % 检查是否陷入局部最优
               if gen > 10 && length(diversity_history) >= 10
                   recent_diversity = mean(diversity_history(end-9:end));
                   if recent_diversity < 0.1  % 多样性过低
                       stagnation_counter = stagnation_counter + 1;
                   else
                       stagnation_counter = 0;
                   end
               end
               
               % 随机重启机制（在 Population 上操作，并保持同步）
               if stagnation_counter >= 5
                   % 随机重启：保留部分优秀解，重新生成部分解
                   [~, sorted_idx] = sort(sum(objvs, 2));
                   keep_count = round(popSize * 0.3);  % 保留30%的优秀解
                   new_count = popSize - keep_count;
                   
                   % 保留优秀解（直接在 Population 上截取）
                   keep_idx   = sorted_idx(1:keep_count);
                   Population = Population(keep_idx);
                   
                   % 生成新的随机解，并加入 Population
                   for i = 1:new_count
                       new_sol = generate_random_solution(varDim);
                       new_sol = repair_route(new_sol);
                       new_obj = Problem.Evaluation(new_sol);
                       Population = [Population, new_obj];
                   end
                   
                   stagnation_counter = 0;
                   diversity_history = [];
                   % 重启后直接进入下一代
                   continue;
               end
               
               % ========== 每代生成 N 个子代 ==========
               % 为每个父代个体生成一个子代（类似 GOCEAfororder.m 的逻辑）
               auxPop = zeros(popSize, varDim);  % 存储所有子代的决策变量
               auxVals = cell(popSize, 1);       % 存储所有子代的 SOLUTION 对象
               
               for i = 1:popSize
                   currentSol = pop(i,:);
                   
                   % 随机选择另一个父代（可以是当前个体自己，也可以是其他个体）
                   idx = randsample(popSize, 1);
                   parent2 = pop(idx,:);
                   
                   % 生成子代（交叉 + 变异 + 修复）
                   trialSol = pop(1,:);  % 临时初始值
                   while any(all(trialSol == pop, 2))
                       child = genetic_crossover_0831(currentSol, parent2);
                       child = bianyi_0831(child);
                       child = repair_route(child);  % 统一修复，确保编码合法
                       trialSol = child;
                   end
                   
                   % 评价子代
                   trialVal = Problem.Evaluation(trialSol);
                   auxPop(i,:) = trialSol;
                   auxVals{i} = trialVal;
               end
               
               % ========== 环境选择：从 2N 个个体中选择 N 个 ==========
               % 合并父代和子代
               PopulationNew = [Population, auxVals{:}];  % 2N 个个体
               % 提取子代的目标值
               auxObjvs = zeros(popSize, objDim);
               for i = 1:popSize
                   auxObjvs(i,:) = auxVals{i}.objs;
               end
               objvsNew = [objvs; auxObjvs];  % 2N 个目标值
               
               % 非支配排序
               [rk,~] = NDSort(objvsNew, inf);
               refPoint = max(objvsNew, [], 1);
               refPoint = 1.2 * refPoint;
               
               % 选择策略：从 2N 个个体中选择 N 个
               selected = false(size(objvsNew,1), 1);
               selected_count = 0;
               
               % 按前沿等级依次选择，直到选够 N 个
               for front = 1:max(rk)
                   front_mask = (rk == front);
                   front_size = sum(front_mask);
                   
                   if selected_count + front_size <= popSize
                       % 整个前沿都能选入
                       selected(front_mask) = true;
                       selected_count = selected_count + front_size;
                   else
                       % 只能选部分：用面积贡献度选择
                       need = popSize - selected_count;
                       FiObjvs = objvsNew(front_mask, 1:objDim);
                       FiPop = PopulationNew(front_mask);
                       
                       [frontObjvs, IX] = sortrows(FiObjvs, 1);
                       fitV = zeros(front_size, 1);
                       
                       if front_size == 1
                           fitV(IX(1)) = 1;  % 只有一个解，直接选
                       else
                           fitV(IX(1)) = (frontObjvs(2,1) - frontObjvs(1,1)) * (refPoint(1,2) - frontObjvs(1,2));
                           if front_size > 2
                               fitV(IX(2:front_size-1)) = (frontObjvs(3:front_size,1) - frontObjvs(2:front_size-1,1)) .* ...
                                                           (frontObjvs(1:front_size-2,2) - frontObjvs(2:front_size-1,2));
                           end
                           fitV(IX(front_size)) = (refPoint(1,1) - frontObjvs(front_size,1)) * ...
                                                   (frontObjvs(front_size-1,2) - frontObjvs(front_size,2));
                       end
                       
                       % 选择面积贡献最大的 need 个
                       [~, sort_idx] = sort(fitV, 'descend');
                       selected_idx = find(front_mask);
                       selected(selected_idx(sort_idx(1:need))) = true;
                       selected_count = popSize;
                   end
                   
                   if selected_count >= popSize
                       break;
                   end
               end
               
               % 更新 Population
               Population = PopulationNew(selected);
           end
        end
    end
end

% 计算种群多样性的辅助函数
function diversity = calculate_diversity(pop, objvs)
    % 基于目标空间的多样性计算
    if size(objvs, 1) < 2
        diversity = 1;
        return;
    end
    
    % 计算目标空间的中心点
    center = mean(objvs, 1);
    
    % 计算每个解到中心的距离
    distances = sqrt(sum((objvs - center).^2, 2));
    
    % 多样性指标：距离的标准差
    diversity = std(distances);
    
    % 归一化到[0,1]范围
    if diversity > 0
        diversity = min(1, diversity / max(distances));
    end
end

% 生成随机解的辅助函数
function sol = generate_random_solution(varDim)
    % 生成1×12的随机解
    sol = zeros(1, varDim);
    
    % 随机选择10个客户编号（1-10）
    customers = randperm(10, 10);
    
    % 随机选择2个0的位置（2-11）
    zero_positions = randperm(10, 2) + 1;
    
    % 组装解
    sol(zero_positions) = 0;
    non_zero_pos = setdiff(1:varDim, zero_positions);
    sol(non_zero_pos) = customers;
end

function d=NumberOfDominatingPoints(P)
[PSize,objDim]=size(P);
d=zeros(PSize,1);
% Calculate number of points from P that dominate each solution in S
for i=1:PSize
    sol=P(i,1:objDim);
    repObjv=sol(ones(PSize,1),1:objDim);                                   % copy current individual
    pos=(P<=repObjv);
    pos=sum(pos,2);
    d(i)=sum(pos==objDim)-1;
end
end

