classdef simplePRNSGAIII < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm III

%------------------------------- Reference --------------------------------
% K. Deb and H. Jain, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part I:
% Solving problems with box constraints, IEEE Transactions on Evolutionary
% Computation, 2014, 18(4): 577-601.
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
            %% Generate the reference points and random population
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
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
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
                
                PopDec=Population.decs;
                varDim=size(PopDec,2);
                for j=1:size(PopDec,1)
                [validPoints]=checknode_singlepath_close(PopDec(j,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei);%对当前解包含的航迹点进行碰撞检测      
               % [PopDecchange,validPoints,bj]=checknode_singlepath_close_newbj(PopDec(j,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);%对当前解包含的航迹点进行碰撞检测
               % Population(j).dec=PopDecchange;
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