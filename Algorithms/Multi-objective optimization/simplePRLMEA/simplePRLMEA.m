classdef simplePRLMEA < ALGORITHM
% <multi/many> <real/integer> <large>
% Evolutionary algorithm for large-scale many-objective optimization
% nSel ---  5 --- Number of selected solutions for decision variable clustering
% nPer --- 50 --- Number of perturbations on each solution for decision variable clustering
% nCor ---  5 --- Number of selected solutions for decision variable interaction analysis
% type ---  1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% X. Zhang, Y. Tian, R. Cheng, and Y. Jin, A decision variable clustering
% based evolutionary algorithm for large-scale many-objective optimization,
% IEEE Transactions on Evolutionary Computation, 2018, 22(1): 97-112.
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
            %% Parameter setting
            [nSel,nPer,nCor,type] = Algorithm.ParameterSet(5,50,5,1);

            %% Generate random population
            Population = Problem.Initialization();
            
            % 小场景环境数据（与 planner_simple_* 配套）
            filename   = 'city_environment_data_new14_simple.xlsx';
            buildings  = xlsread(filename,'Buildings');
            numBuildings = size(buildings,1);
            cylinders  = xlsread(filename,'Cylinders');
            numCylinders = size(cylinders,1);
            pyramids   = xlsread(filename,'Pyramids');
            numPyramids = size(pyramids,1);
            spheres    = xlsread(filename,'Spheres');
            numSpheres = size(spheres,1);
            jinfei     = xlsread(filename,'jinfei');
            numjinfei  = size(jinfei,1);

            %% Detect the group of each distance variable
            [PV,DV] = VariableClustering(Problem,Population,nSel,nPer);
            DVSet   = CorrelationAnalysis(Problem,Population,DV,nCor);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Convergence optimization
                for i = 1 : 10
                    drawnow('limitrate');
                    Population = ConvergenceOptimization(Problem,Population,DVSet,type);
                end
                % Distribution optimization
                for i = 1 : Problem.M
                    drawnow();
                    Population = DistributionOptimization(Problem,Population,PV);
                end
                
                % 碰撞检测，标记 validtrait（1=不碰撞，0=碰撞）
                PopDec = Population.decs;
                varDim = size(PopDec,2);
                for j = 1:size(PopDec,1)
                    validPoints = checknode_singlepath_close(PopDec(j,:),varDim, ...
                        buildings,numBuildings,cylinders,numCylinders, ...
                        spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei);
                    if any(~validPoints)
                        Population(j).validtrait = 0;
                    else
                        Population(j).validtrait = 1;
                    end
                end
            end
        end
    end
end

function Population = ConvergenceOptimization(Problem,Population,DVSet,type)
    [N,D] = size(Population.decs);
    Con   = calCon(Population.objs);
    for i = 1 : length(DVSet)
        for j = 1 : length(DVSet{i})
            % Select parents
            MatingPool = TournamentSelection(2,2*N,Con);
            % Generate offsprings
            OffDec = Population.decs;
            if type == 1
                NewDec = OperatorGAhalf(Problem,Population(MatingPool).decs,{1,20,D/length(DVSet{i})/2,20});
            elseif type == 2
                NewDec = OperatorDE(Problem,Population.decs,Population(MatingPool(1:end/2)).decs,Population(MatingPool(end/2+1:end)).decs,{1,0.5,D/length(DVSet{i})/2,20});
            end
            OffDec(:,DVSet{i}) = NewDec(:,DVSet{i});
            Offspring          = Problem.Evaluation(OffDec);
            % Update each solution
            allCon  = calCon([Population.objs;Offspring.objs]);
            Con     = allCon(1:N);
            newCon  = allCon(N+1:end);
            updated = Con > newCon;
            Population(updated) = Offspring(updated);
            Con(updated)        = newCon(updated);
        end
    end
end

function Population = DistributionOptimization(Problem,Population,PV)
% Distribution optimization

    N            = length(Population);
    OffDec       = Population(TournamentSelection(2,N,calCon(Population.objs))).decs;
    NewDec       = OperatorGA(Problem,Population(randi(N,1,N)).decs);
    OffDec(:,PV) = NewDec(:,PV);
    Offspring    = Problem.Evaluation(OffDec);
    Population   = EnvironmentalSelection([Population,Offspring],N);
end

function Con = calCon(PopObj)
% Calculate the convergence of each solution

    FrontNo = NDSort(PopObj,inf);
    Con     = sum(PopObj,2);
    Con     = FrontNo'*(max(Con)-min(Con)) + Con;
end