classdef SMSEMOA_G < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% SMS-EMOA (OperatorAMO variant)
%
% 参数说明：
% pc --- 0.5  --- 交叉概率（交叉/变异二选一）
% er --- 0.25 --- explorationRate（变异算子中 T4 大步长概率基值）
%
%------------------------------- Reference --------------------------------
% M. Emmerich, N. Beume, and B. Naujoks, An EMO algorithm using the
% hypervolume measure as selection criterion, Proceedings of the
% International Conference on Evolutionary Multi-Criterion Optimization,
% 2005, 62-76.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameters
            [pc,er] = Algorithm.ParameterSet(0.5,0.25);

            %% Generate random population
            Population = Problem.Initialization();
            FrontNo    = NDSort(Population.objs,inf);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                for i = 1 : Problem.N
                    drawnow('limitrate');
                    idx = randperm(length(Population),2);
                    Parent1 = Population(idx(1));
                    Parent2 = Population(idx(2));
                    Offspring = OperatorAMO(Problem,Parent1,'P2',Parent2,'pc',pc,'er',er);

                    % Decision-space deduplication (skip clones)
                    try
                        if ismember(Offspring.dec,Population.decs,'rows')
                            continue;
                        end
                    catch
                    end

                    [Population,FrontNo] = Reduce([Population,Offspring],FrontNo);
                end
            end
        end
    end
end

