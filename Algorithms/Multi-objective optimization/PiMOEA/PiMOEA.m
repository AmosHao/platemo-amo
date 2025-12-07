classdef PiMOEA < ALGORITHM
% <many> <real/binary/permutation> <constrained/none>
% Pivot-based Evolutionary algorithm
% 
%------------------------------- Reference --------------------------------
% Palakonda, Vikas, Jae-Mo Kang, and Heechul Jung. "An Adaptive Neighborhood based 
% Evolutionary Algorithm with Pivot-Solution based Selection for Multi-and 
% Many-Objective Optimization." Information Sciences (2022).
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting           

            %% Generate random population
            Population = Problem.Initialization();
            FrontNo    = NDSort(Population.objs,Population.cons,inf);
            PivotSol   = zeros(1,Problem.N);     % Set of pivot solutions
            r          = -ones(1,2*Problem.N);	% Ratio of size of neighorhood
            t          = -ones(1,2*Problem.N);	% Ratio of pivot solutions

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs,FrontNo,PivotSol);
                %Offspring  = OperatorGA(Population(MatingPool));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = [Population,Offspring];
                [FrontNo,MaxFNo]              = NDSort(Population.objs,Population.cons,Problem.N);
                [PivotSol,Distance,r,t]       = F_PivotSol(Population.objs,FrontNo,MaxFNo,r,t);
                [Population,FrontNo,PivotSol] = EnvironmentalSelection(Population,FrontNo,MaxFNo,PivotSol,Distance,Problem.N);      
            end
        end
    end
end