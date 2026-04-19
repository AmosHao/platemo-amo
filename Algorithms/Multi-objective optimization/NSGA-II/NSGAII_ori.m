classdef NSGAII_ori < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none>
% NSGA-II baseline with RepairOrderEncoding (dec -> repair -> evaluation)
%
%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAII(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);

                ParentsDec = Population(MatingPool).decs;
                OffDec     = OperatorGA(Problem,ParentsDec);
                OffDec     = RepairOrderEncoding(Problem,OffDec);
                Offspring  = Problem.Evaluation(OffDec);

                [Population,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAII([Population,Offspring],Problem.N);
            end
        end
    end
end

