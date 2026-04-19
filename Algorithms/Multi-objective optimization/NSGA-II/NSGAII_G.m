classdef NSGAII_G < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none>
% NSGA-II (OperatorAMO variant)
%
% 参数说明：
% pc --- 0.5  --- 交叉概率（交叉/变异二选一）
% er --- 0.25 --- explorationRate（变异算子中 T4 大步长概率基值）
%
%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameters
            [pc,er] = Algorithm.ParameterSet(0.5,0.25);

            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAII(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Parents    = Population(MatingPool);
                Offspring  = OperatorAMO(Problem,Parents,'pc',pc,'er',er);

                Union = DeduplicateByDec([Population,Offspring]);
                if length(Union) >= Problem.N
                    [Population,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAII(Union,Problem.N);
                else
                    [Population,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAII([Population,Offspring],Problem.N);
                end
            end
        end
    end
end

