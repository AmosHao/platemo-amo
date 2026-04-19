classdef PeEA_ori < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% PeEA baseline with RepairOrderEncoding (dec -> repair -> evaluation)
%
%------------------------------- Reference --------------------------------
% L. Li, G. G. Yen, A. Sahoo, L. Chang, and T. Gu, On the estimation of
% pareto front and dimensional similarity in many-objective evolutionary
% algorithm, Information Sciences, 2021, 563: 375-400.
%--------------------------------------------------------------------------

	methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs);

                ParentsDec = Population(MatingPool).decs;
                OffDec     = OperatorGA(Problem,ParentsDec);
                OffDec     = RepairOrderEncoding(Problem,OffDec);
                Offspring  = Problem.Evaluation(OffDec);

                Population = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end

