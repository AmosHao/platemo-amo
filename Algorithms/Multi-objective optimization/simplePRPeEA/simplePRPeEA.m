classdef simplePRPeEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Pareto front shape estimation based evolutionary algorithm

%------------------------------- Reference --------------------------------
% L. Li, G. G. Yen, A. Sahoo, L. Chang, and T. Gu, On the estimation of
% pareto front and dimensional similarity in many-objective evolutionary
% algorithm, Information Sciences, 2021, 563: 375-400.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Li Li

	methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
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
                MatingPool = MatingSelection(Population.objs);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = EnvironmentalSelection([Population,Offspring],Problem.N); 
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