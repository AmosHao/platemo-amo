classdef simplePRPiMOEA < ALGORITHM
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
                MatingPool = MatingSelection(Population.objs,FrontNo,PivotSol);
                %Offspring  = OperatorGA(Population(MatingPool));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = [Population,Offspring];
                [FrontNo,MaxFNo]              = NDSort(Population.objs,Population.cons,Problem.N);
                [PivotSol,Distance,r,t]       = F_PivotSol(Population.objs,FrontNo,MaxFNo,r,t);
                [Population,FrontNo,PivotSol] = EnvironmentalSelection(Population,FrontNo,MaxFNo,PivotSol,Distance,Problem.N);      
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