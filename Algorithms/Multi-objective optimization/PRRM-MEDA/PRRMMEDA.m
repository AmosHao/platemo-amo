classdef PRRMMEDA < ALGORITHM
% <multi> <real/integer>
% Regularity model-based multiobjective estimation of distribution
% algorithm (Path Repair)

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, and Y. Jin, RM-MEDA: A regularity model-based
% multiobjective estimation of distribution algorithm, IEEE Transactions on
% Evolutionary Computation, 2008, 12(1): 41-63.
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
            %% Generate random population
            Population = Problem.Initialization();

            %% 场景设置
            buildings = load('buildings.mat').buildings;
            numBuildings = size(buildings, 1);
            cylinders = load('cylinders.mat').cylinders;
            numCylinders = size(cylinders, 1);
            pyramids = load('pyramids.mat').pyramids;
            numPyramids = size(pyramids, 1);
            spheres = load('spheres.mat').spheres;
            numSpheres = size(spheres, 1);
            jinfei = load('jinfei.mat').jinfei;
            numjinfei = size(jinfei, 1);
            dots = Problem.dots;
            dot1 = [dots(1,1),dots(1,2),dots(1,3)];
            dot2 = [dots(1,4),dots(1,5),dots(1,6)];

            %% Optimization
            while Algorithm.NotTerminated(Population)
                Offspring  = Operator(Problem,Population);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N);

                %% 路径修复与碰撞检测
                PopDec = Population.decs;
                varDim = size(PopDec,2);
                for j = 1:size(PopDec,1)
                    [PopDecchange,validPoints,bj] = checknode_singlepath_close_newbj(PopDec(j,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);
                    allVals = Problem.Evaluation(PopDecchange);
                    Population(j) = allVals;
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
