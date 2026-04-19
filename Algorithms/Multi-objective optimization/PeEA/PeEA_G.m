classdef PeEA_G < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Pareto front shape estimation based evolutionary algorithm (OperatorAMO variant)
%
% 参数说明：
% pc --- 0.5  --- 交叉概率（交叉/变异二选一）
% er --- 0.25 --- explorationRate（变异算子中 T4 大步长概率基值）
%
%------------------------------- Reference --------------------------------
% L. Li, G. G. Yen, A. Sahoo, L. Chang, and T. Gu, On the estimation of
% pareto front and dimensional similarity in many-objective evolutionary
% algorithm, Information Sciences, 2021, 563: 375-400.
%--------------------------------------------------------------------------

	methods
        function main(Algorithm,Problem)
            %% Parameters
            [pc,er] = Algorithm.ParameterSet(0.5,0.25);

            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs);
                Parents    = Population(MatingPool);
                Offspring  = OperatorAMO(Problem,Parents,'pc',pc,'er',er);
                Union = DeduplicateByDec([Population,Offspring]);
                if length(Union) >= Problem.N
                    Population = EnvironmentalSelection(Union,Problem.N);
                else
                    Population = EnvironmentalSelection([Population,Offspring],Problem.N);
                end
            end
        end
    end
end

