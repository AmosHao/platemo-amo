classdef GOCEAfororder_0722< ALGORITHM
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
               
               idx      = randsample(popSize,2); 
               parents  = pop(idx,:);
               trialSol = pop(1,:);
               while any(all(trialSol == pop, 2))
                   % trialSol=DifferentialEvolutionCrossover(parents,bounds,F,CR); 
                   % [child1, child2]=genetic_crossover_0831(parents(1,:),parents(2,:));
                   % % trialSol=PolynomialMutation(trialSol,bounds,pm);
                   % childs=[child1;child2];
                   % index=randi([1,2]);
                   % child=childs(index,:);
                  child = genetic_crossover_0831(parents(1,:), parents(2,:));   % 直接返回一条 12 维
                  child = bianyi_0831(child);
                  % 统一修复，确保编码合法（10个客户+2个0）
                  child = repair_route(child);
                  trialSol = child;
                end
               % 评价新解，并与当前 Population 合并
               trialVal      = Problem.Evaluation(trialSol);
               PopulationNew = [Population, trialVal];
               objvsNew      = [objvs; trialVal.objs];
               
               % 基于合并后的种群进行非支配排序和删选
               [rk,~]   = NDSort(objvsNew,inf);
               refPoint = max(objvsNew,[],1); 
               refPoint = 1.2*refPoint;

               if max(rk) ~= 1
                   % 存在多个前沿：删掉“被支配最多”的个体
                   d = NumberOfDominatingPoints(objvsNew); 
                   [~,loc] = max(d);
                   keep = true(size(objvsNew,1),1);
                   keep(loc) = false;
                   Population = PopulationNew(keep);
               else
                   % 全部在同一前沿：用面积近似贡献度来删一个
                   FiMask  = (rk == max(rk));
                   FiSize  = sum(FiMask);
                   FiObjvs = objvsNew(FiMask,1:objDim);
                   FiPop   = PopulationNew(FiMask);
                   
                   frontObjvs = FiObjvs;
                   fitV       = zeros(FiSize,1);
                   [frontObjvs,IX] = sortrows(frontObjvs,1);
                   
                   fitV(IX(1))          = (frontObjvs(2,1)-frontObjvs(1,1))   .* (refPoint(1,2)-frontObjvs(1,2));
                   fitV(IX(2:FiSize-1)) = (frontObjvs(3:FiSize,1)-frontObjvs(2:FiSize-1,1)) .* (frontObjvs(1:FiSize-2,2)-frontObjvs(2:FiSize-1,2));
                   fitV(IX(FiSize))     = (refPoint(1,1)-frontObjvs(FiSize,1)) .* (frontObjvs(FiSize-1,2)-frontObjvs(FiSize,2));
                   
                   [~,locFi] = min(fitV);
                   FiObjvs(locFi,:) = [];
                   FiPop(locFi)     = [];
                   
                   % 重新组合各前沿，得到新的 Population
                   Population = [PopulationNew(~FiMask), FiPop];
               end
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