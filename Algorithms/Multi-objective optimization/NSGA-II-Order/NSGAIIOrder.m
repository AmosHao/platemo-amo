classdef NSGAIIOrder < ALGORITHM
% <multi> <real/integer> <constrained/none>
% NSGA-II for Order Assignment Problem
% 专门适配planner_order_v3问题的NSGA-II变体
% 使用自定义的交叉和变异操作符处理带0分隔符的序列

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAIIOrder(Population,Problem.N);

            %% 多样性监控和自适应参数
            gen = 0;
            diversity_history = [];
            stagnation_counter = 0;
            maxgen = 400;  % 最大迭代次数
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen = gen + 1;
                
                % 计算当前多样性（基于目标空间）
                if gen > 1
                    objVals = Population.objs;
                    % 计算目标值的标准差作为多样性指标
                    if size(objVals, 1) > 1
                        diversity = mean(std(objVals, 0, 1)) / (mean(mean(objVals)) + 1e-10);
                        diversity_history = [diversity_history, diversity];
                        
                        % 检测停滞
                        if gen > 20 && length(diversity_history) >= 20
                            recent_diversity = mean(diversity_history(end-19:end));
                            if recent_diversity < 0.01  % 多样性过低
                                stagnation_counter = stagnation_counter + 1;
                            else
                                stagnation_counter = max(0, stagnation_counter - 1);
                            end
                        end
                    end
                end
                
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorOrder(Problem,Population(MatingPool),gen,maxgen,stagnation_counter);
                [Population,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAIIOrder([Population,Offspring],Problem.N);
                
                % 如果停滞，增加随机重启
                if stagnation_counter >= 10 && mod(gen, 10) == 0
                    % 随机替换10%的种群
                    numReplace = max(1, floor(Problem.N * 0.1));
                    replaceIdx = randperm(Problem.N, numReplace);
                    newPop = Problem.Initialization(numReplace);
                    for i = 1:numReplace
                        Population(replaceIdx(i)) = newPop(i);
                    end
                    [~,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAIIOrder(Population,Problem.N);
                    stagnation_counter = 0;  % 重置计数器
                end
            end
        end
    end
end

function Offspring = OperatorOrder(Problem,Parent,gen,maxgen,stagnation_counter)
% 自定义操作符：处理带0分隔符的序列
% 使用A-amos算法中的交叉和变异操作
% 增加自适应参数和多样性增强机制

    if nargin < 3
        gen = 1;
        maxgen = 400;
        stagnation_counter = 0;
    end
    
    if isa(Parent(1),'SOLUTION')
        evaluated = true;
        ParentDec = Parent.decs;
    else
        evaluated = false;
        ParentDec = Parent;
    end
    
    [N,D] = size(ParentDec);
    Parent1 = ParentDec(1:floor(end/2),:);
    Parent2 = ParentDec(floor(end/2)+1:floor(end/2)*2,:);
    Offspring = zeros(2*size(Parent1,1),D);
    
    % 获取问题参数
    n = Problem.n;  % 客户点数量（商家数量）
    m = Problem.m;  % 无人机数量
    
    % 自适应交叉概率：早期高交叉，后期降低
    proC = 0.9 - 0.2 * (gen / maxgen);
    
    % 自适应变异概率：根据停滞情况动态调整
    base_pm = 0.15;
    if stagnation_counter > 5
        % 如果停滞，大幅增加变异概率
        proM = min(0.5, base_pm * (1 + stagnation_counter * 0.1));
    else
        % 早期高变异探索，后期降低
        proM = base_pm * (1.5 - 0.5 * gen / maxgen);
    end
    
    % 添加A-amos算法文件夹到路径（如果函数不在当前路径）
    aamosPath = fileparts(mfilename('fullpath'));
    aamosPath = fullfile(fileparts(aamosPath), 'A-amos');
    if exist(aamosPath, 'dir')
        addpath(aamosPath);
    end
    
    % 对每一对父代进行交叉和变异
    for i = 1:size(Parent1,1)
        % 确保维度正确
        p1 = Parent1(i,:);
        p2 = Parent2(i,:);
        if length(p1) ~= D
            if length(p1) < D
                p1 = [p1, zeros(1, D-length(p1))];
            else
                p1 = p1(1:D);
            end
        end
        if length(p2) ~= D
            if length(p2) < D
                p2 = [p2, zeros(1, D-length(p2))];
            else
                p2 = p2(1:D);
            end
        end
        
        % 交叉操作
        if rand < proC
            % 使用A-amos的段级交叉
            try
                child = SegmentBasedCrossover(p1, p2, m);
            catch ME
                % 如果交叉失败，随机选择一个父代
                child = p1;
            end
        else
            child = p1;
        end
        
        % 确保维度
        if length(child) ~= D
            if length(child) < D
                child = [child, zeros(1, D-length(child))];
            else
                child = child(1:D);
            end
        end
        
        % 变异操作（增强：多次变异以增加多样性）
        if rand < proM
            % 使用A-amos的离散变异
            try
                child = DiscreteMutation(child, m);
                % 如果停滞，增加额外的变异
                if stagnation_counter > 5 && rand < 0.3
                    child = DiscreteMutation(child, m);
                end
            catch ME
                % 如果变异失败，保持原样
            end
        end
        
        % 增加随机扰动（如果停滞）
        if stagnation_counter > 10 && rand < 0.2
            % 随机交换两个非0位置
            nonZeroIdx = find(child ~= 0);
            if length(nonZeroIdx) >= 2
                swapIdx = randperm(length(nonZeroIdx), 2);
                temp = child(nonZeroIdx(swapIdx(1)));
                child(nonZeroIdx(swapIdx(1))) = child(nonZeroIdx(swapIdx(2)));
                child(nonZeroIdx(swapIdx(2))) = temp;
            end
        end
        
        % 确保维度
        if length(child) ~= D
            if length(child) < D
                child = [child, zeros(1, D-length(child))];
            else
                child = child(1:D);
            end
        end
        
        % 解的完整性修复（关键：确保包含所有必需的点）
        try
            child = EnsureSolutionCompleteness(child, n, m, D);
        catch ME
            % 如果修复失败，尝试重新生成
            try
                child = Problem.Initialization(1).decs(1,:);
            catch
            end
        end
        
        % 顺序约束修复
        try
            child = OrderConstraintRepair(child, n, m);
        catch ME
            % 如果修复失败，保持原样
        end
        
        % 最终确保维度
        if length(child) ~= D
            if length(child) < D
                child = [child, zeros(1, D-length(child))];
            else
                child = child(1:D);
            end
        end
        
        Offspring(i,:) = child;
        
        % 第二个子代（从另一个父代开始）
        if rand < proC
            try
                child2 = SegmentBasedCrossover(p2, p1, m);
            catch ME
                child2 = p2;
            end
        else
            child2 = p2;
        end
        
        % 确保维度
        if length(child2) ~= D
            if length(child2) < D
                child2 = [child2, zeros(1, D-length(child2))];
            else
                child2 = child2(1:D);
            end
        end
        
        if rand < proM
            try
                child2 = DiscreteMutation(child2, m);
                % 如果停滞，增加额外的变异
                if stagnation_counter > 5 && rand < 0.3
                    child2 = DiscreteMutation(child2, m);
                end
            catch ME
            end
        end
        
        % 增加随机扰动（如果停滞）
        if stagnation_counter > 10 && rand < 0.2
            % 随机交换两个非0位置
            nonZeroIdx = find(child2 ~= 0);
            if length(nonZeroIdx) >= 2
                swapIdx = randperm(length(nonZeroIdx), 2);
                temp = child2(nonZeroIdx(swapIdx(1)));
                child2(nonZeroIdx(swapIdx(1))) = child2(nonZeroIdx(swapIdx(2)));
                child2(nonZeroIdx(swapIdx(2))) = temp;
            end
        end
        
        % 确保维度
        if length(child2) ~= D
            if length(child2) < D
                child2 = [child2, zeros(1, D-length(child2))];
            else
                child2 = child2(1:D);
            end
        end
        
        % 解的完整性修复（关键：确保包含所有必需的点）
        try
            child2 = EnsureSolutionCompleteness(child2, n, m, D);
        catch ME
            % 如果修复失败，尝试重新生成
            try
                child2 = Problem.Initialization(1).decs(1,:);
            catch
            end
        end
        
        try
            child2 = OrderConstraintRepair(child2, n, m);
        catch ME
        end
        
        % 最终确保维度
        if length(child2) ~= D
            if length(child2) < D
                child2 = [child2, zeros(1, D-length(child2))];
            else
                child2 = child2(1:D);
            end
        end
        
        Offspring(i+size(Parent1,1),:) = child2;
    end
    
    if evaluated
        Offspring = Problem.Evaluation(Offspring);
    end
end
