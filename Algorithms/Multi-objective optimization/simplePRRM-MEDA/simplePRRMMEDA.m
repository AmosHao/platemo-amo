classdef simplePRRMMEDA < ALGORITHM
% <multi> <real/integer>
% Regularity model-based multiobjective estimation of distribution
% algorithm

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

            %% Optimization
            while Algorithm.NotTerminated(Population)
                Offspring  = Operator(Problem,Population);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N);
                
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