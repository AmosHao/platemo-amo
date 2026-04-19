classdef SMSEMOA_ori < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% SMS-EMOA baseline with RepairOrderEncoding (dec -> repair -> evaluation)
%
%------------------------------- Reference --------------------------------
% M. Emmerich, N. Beume, and B. Naujoks, An EMO algorithm using the
% hypervolume measure as selection criterion, Proceedings of the
% International Conference on Evolutionary Multi-Criterion Optimization,
% 2005, 62-76.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            FrontNo    = NDSort(Population.objs,inf);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                for i = 1 : Problem.N
                    drawnow('limitrate');
                    idx = randperm(length(Population),2);

                    ParentsDec = Population(idx).decs;
                    OffDec     = OperatorGAhalf(Problem,ParentsDec);   % 2 parents -> 1 offspring dec
                    OffDec     = RepairOrderEncoding(Problem,OffDec);
                    Offspring  = Problem.Evaluation(OffDec);

                    [Population,FrontNo] = Reduce([Population,Offspring],FrontNo);
                end
            end
        end
    end
end

