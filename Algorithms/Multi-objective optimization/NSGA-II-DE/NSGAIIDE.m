classdef NSGAIIDE < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm II

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

        function main(Algorithm,Problem)%main函数
  

            %% Generate random population初 始化种群
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            %bounds = [Problem.lower,Problem.upper];
            %% Optimization
            while Algorithm.NotTerminated(Population)%调用algorithm中的notterminated
                MatingPool = TournamentSelection(3,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorDE(Problem, Population(MatingPool(1)), Population(MatingPool(2)), Population(MatingPool(3)));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end

      
       % function [objvs, pop] = main(~, prob)
       %      maxGens=300;
       %      F= 0.3;
       %      CR = 0.2;
       %      D = prob.D;            % The number of decision variables
       %      M = prob.M;            % The number of objectives 
       %      if M == 2
       %          popSize = 100;       % Population size
       %      else
       %          popSize = 105;       % Population size
       %      end
       % 
       %      % Mutation probability
       %      Pm = 1.0 / D;
       % 
       %      % Set the lower and upper boundaries
       %      bounds = [prob.lower;prob.upper];
       % 
       %      % Initialize population and objective values
       %      [pop,objvs] = InitializationDE(prob,popSize);
       %      objDim = size(objvs, 2);
       %      % Set initial offspring population
       %      offsPop = pop;
       %      offsObjvs = objvs;
       % 
       %      % Preallocation for parents
       %      parents = inf(3, D);             
       % 
       %      for gen = 1:maxGens
       %          for i = 1:popSize
       %              parents(3, :) = pop(i, 1:D);
       %              parents(1:2, :) = BinaryTournamentSelection(2, pop, objvs); % Select two different parents
       % 
       %              offs = DifferentialEvolutionCrossover(parents, bounds, F, CR);
       %              offs = PolynomialMutation(offs(1, 1:D), bounds, Pm);
       % 
       %              offsPop(i, 1:D) = offs;
       %              offsObjvs(i, 1:objDim) = SolutionEvaluationDE(prob, offs);
       %          end
       %      end
       % end
 
       
    end
 
end