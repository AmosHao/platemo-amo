classdef NSGAIII_G < ALGORITHM
% <multi/many> <real/integer>
% NSGA-III-G: 仅嵌入 A-AMOS 的结构保持交叉/变异算子（OperatorAMO），交配选择仍用原锦标赛
%
% 参数说明：
% pc --- 0.5  --- 交叉概率（交叉/变异二选一）
% er --- 0.25 --- explorationRate（变异算子中 T4 大步长概率基值）

    methods
        function main(Algorithm,Problem)
            %% Parameters
            [pc,er] = Algorithm.ParameterSet(0.5,0.25);

            %% Generate the reference points and random population
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % 为 OperatorAMO 提供 N 父代（每个父代生成 1 个子代；交叉时在该集合内随机选另一个父代）
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Parents    = Population(MatingPool);
                Offspring  = OperatorAMO(Problem,Parents,'pc',pc,'er',er);

                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);

                Union = DeduplicateByDec([Population,Offspring]);
                if length(Union) >= Problem.N
                    Population = EnvironmentalSelection(Union,Problem.N,Z,Zmin);
                else
                    Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
                end
            end
        end
    end
end

