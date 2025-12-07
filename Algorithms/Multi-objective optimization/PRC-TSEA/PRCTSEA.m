classdef PRCTSEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained>
% Constrained two-stage evolutionary algorithm

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, H. Zhen, S. Li, L. Wang, and Z. Liao, A simple
% two-stage evolutionary algorithm for constrained multi-objective
% optimization, Knowledge-Based Systems, 2021, 228: 107263.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    methods
        function main(Algorithm,Problem)
            %% Generate the sampling points and random population
            Population = Problem.Initialization();
            W = UniformPoint(Problem.N,Problem.M);
            [ARMOEA_Archive,RefPoint,Range] = UpdateRefPoint(Population.objs,W,[]);
            CV = sum(max(0,Population.cons),2);
            Archive = Population(CV==0);
            stage_conver = 0;
             %场景设置
        buildings = load('buildings.mat').buildings;
        numBuildings = size(buildings, 1);
        % 读取圆柱体数据
       cylinders = load('cylinders.mat').cylinders;
        numCylinders = size(cylinders, 1);
        % 读取四棱锥数据
        pyramids = load('pyramids.mat').pyramids;
        numPyramids = size(pyramids, 1);
        % 读取球体数据
        spheres = load('spheres.mat').spheres;
        numSpheres = size(spheres, 1);
        % 读取禁飞区数据
        jinfei = load('jinfei.mat').jinfei;
        numjinfei = size(jinfei, 1);
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Problem.FE<0.5*Problem.maxFE
                    % evolve population to PF by ARMOEA
                    MatingPool = MatingSelection1(Population,RefPoint,Range);
                    Offspring = OperatorGA(Problem,Population(MatingPool));
                    [ARMOEA_Archive,RefPoint,Range] = UpdateRefPoint([ARMOEA_Archive;Offspring.objs],W,Range);
                    Archive = UpdateArchive(Archive,[Population,Offspring],Problem.N);
                    [Population,Range] = EnvironmentalSelection1([Population,Offspring],RefPoint,Range,Problem.N);
                else
                    if stage_conver==0
                        % exchange archive and population
                        temp = Population;
                        Population = Archive;
                        Archive = temp;
                        stage_conver = 1;
                    end
                    % evolve population to CPF by modified SPEA2
                    MatingPool = MatingSelection2(Population,Archive,Problem.N);
                    Offspring = OperatorGA(Problem,MatingPool);
                    Population = EnvironmentalSelection2([Population,Offspring],Problem.N);

                end
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