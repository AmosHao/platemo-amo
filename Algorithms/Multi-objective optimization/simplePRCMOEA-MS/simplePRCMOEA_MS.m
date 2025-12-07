classdef simplePRCMOEA_MS < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained multiobjective evolutionary algorithm with multiple stages
% type   ---   1 --- Type of operator (1. GA 2. DE)
% lambda --- 0.5 --- Parameter for determining the current stage

%------------------------------- Reference --------------------------------
% Y. Tian, Y. Zhang, Y. Su, X. Zhang, K. C. Tan, and Y. Jin, Balancing
% objective optimization and constraint satisfaction in constrained
% evolutionary multi-objective optimization, IEEE Transactions on
% Cybernetics, 2022, 52(9): 9559-9572.
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
        	[type,lambda] = Algorithm.ParameterSet(1,0.5);

           %% Generate random population
            Population = Problem.Initialization();
            Fitness    = CalFitness([CalSDE(Population.objs)',CalCV(Population.cons)]);
% 读取cubes数据
             % 文件名
            filename = 'city_environment_data_new14_simple.xlsx';
            % 读取长方体数据
            buildings = xlsread(filename, 'Buildings');
            numBuildings = size(buildings, 1);
            % 读取圆柱体数据
            cylinders = xlsread(filename, 'Cylinders');
            numCylinders = size(cylinders, 1);
            % 读取四棱锥数据
            pyramids = xlsread(filename, 'Pyramids');
            numPyramids = size(pyramids, 1);
            % 读取球体数据
            spheres = xlsread(filename, 'Spheres');
            numSpheres = size(spheres, 1);
            % 读取禁飞区数据
            jinfei = xlsread(filename, 'jinfei');
            numjinfei = size(jinfei, 1);
            %% Optimization
            while Algorithm.NotTerminated(Population)
                
                if type == 1
                    MatingPool = TournamentSelection(2,Problem.N,Fitness);  
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                elseif type == 2
                    Mat1 = TournamentSelection(2,Problem.N,Fitness);
                    Mat2 = TournamentSelection(2,Problem.N,Fitness);
                    Offspring = OperatorDE(Problem,Population,Population(Mat1),Population(Mat2));
                end
                Q  = [Population,Offspring];
                CV = CalCV(Q.cons);
                if mean(CV<=0) > lambda && Problem.FE >= 0.1*Problem.maxFE
                    Fitness = CalFitness(Q.objs,CV);
                else
                    Fitness = CalFitness([CalSDE(Q.objs)',CV]);
                end
                [Population,Fitness] = EnvironmentalSelection(Fitness,Q,Problem.N);
            PopDec=Population.decs;
                varDim=size(PopDec,2);
                for j=1:size(PopDec,1)
                [validPoints]=checknode_singlepath_close(PopDec(j,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei);%对当前解包含的航迹点进行碰撞检测      
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

function CV = CalCV(CV_Original)
    CV_Original = max(CV_Original,0);
    CV = CV_Original./max(CV_Original,[],1);
    CV(:,isnan(CV(1,:))) = 0;
    CV = mean(CV,2);
end