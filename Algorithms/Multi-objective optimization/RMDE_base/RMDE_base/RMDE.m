classdef RMDE< ALGORITHM
% <multi> <real/integer>
% GOCEA
% Kmax ---   15 --- The K
% BETA    --- 0.5 --- The B
% F ---   0.5 --- The F
% CR    --- 1 --- The CR
    methods
        function main(Algorithm, Problem)
% <algorithm> <R>

% K --- 10 --- ¾ÛÀàÊýÄ¿

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, and Y. Jin, RM-MEDA: A regularity model-based
% multiobjective estimation of distribution algorithm, IEEE Transactions on
% Evolutionary Computation, 2008, 12(1): 41-63.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% ???°è?¾ç½®
    [K] = Algorithm.ParameterSet(10);

    %% Generate random population
    Population = Problem.Initialization();
    
    %% Optimization
    
    while Algorithm.NotTerminated(Population)
        Offspring  = Operator(Population,K,Problem);
	Population = EnvironmentalSelection([Population,Offspring],Problem.N);
    end
        end
    end
end