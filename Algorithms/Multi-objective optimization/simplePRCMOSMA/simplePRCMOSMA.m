classdef simplePRCMOSMA < ALGORITHM
% <multi/many> <real/integer> <constrained>
% Constrained multi-objective evolutionary algorithm with self-organizing map

%------------------------------- Reference --------------------------------
% C. He, M. Li, C. Zhang, H. Chen, P. Zhong, Z. Li, and J. Li, A
% self-organizing map approach for constrained multi-objective optimization
% problems, Complex & Intelligent Systems, 2022, 8: 5355-5375.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Chao He

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [D,tau0,H] = Algorithm.ParameterSet(repmat(ceil(Problem.N.^(1/(Problem.M-1))),1,Problem.M-1),0.7,5);
            D1         = D;
            Problem.N  = prod(D);
            sigma0     = sqrt(sum(D.^2)/(Problem.M-1))/2;

            %% Generate random population
            FP = Problem.Initialization();
            AP = Problem.Initialization();

            %% Initialize the SOM
            % Training set
            S = FP.decs;
            % Weight vector of each neuron
            W = S;
            [LDis,B] = Initialize_SOM(S,D,H);
            % Training set
            S2 = AP.decs;
            % Weight vector of each neuron
            W2 = S2;
            [LDis2,B2] = Initialize_SOM(S2,D1,H);
                         %场景设置
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
            while Algorithm.NotTerminated(FP)
                % Update SOM1
                W = UpdateSOM(S,W,Problem.FE,Problem.maxFE,LDis,sigma0,tau0); 
                % Update SOM2
                W2 = UpdateSOM(S2,W2,Problem.FE,Problem.maxFE,LDis2,sigma0,tau0);

                % Associate each solution with a neuron
                XU  = Associate(FP,W,Problem.N);
                XU2 = Associate(AP,W2,Problem.N); 

                %Construct  matingPool  
                [MatingPool1] = MatingPool(XU,Problem.N,B);
                [MatingPool2] = MatingPool(XU2,Problem.N,B2);

                % Evolution
                A1 = FP.decs;
                Offspring1 = OperatorGA(Problem,[FP(XU),FP(MatingPool1)]);%GA
                A2         = AP.decs;
                Offspring2 = OperatorGA(Problem,[AP(XU2),AP(MatingPool2)]);%GA
                %EnvironmentalSelection
                FP = EnvironmentalSelection([FP,Offspring1,Offspring2],Problem.N,true);
                AP = EnvironmentalSelection([AP,Offspring1,Offspring2],Problem.N,false);
                % Update the training set
                S  = setdiff(FP.decs,A1,'rows');
                S2 = setdiff(AP.decs,A2,'rows');

                PopDec=FP.decs;
                varDim=size(PopDec,2);
                for j=1:size(PopDec,1)
                [validPoints]=checknode_singlepath_close(PopDec(j,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei);%对当前解包含的航迹点进行碰撞检测      
                if any(~validPoints) 
                FP(j).validtrait = 0;
                else   
                FP(j).validtrait = 1;
                end
                end
            end
        end
    end
end