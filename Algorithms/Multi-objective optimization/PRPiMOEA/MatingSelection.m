function MatingPool = MatingSelection(PopObj,FrontNo,PivotSol)
% The mating selection of KnEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the weighted distance of each solution   
    Crowd             = F_AVRank(PopObj);

    %% Binary tournament selection
    MatingPool = TournamentSelection(2,size(PopObj,1),FrontNo,-PivotSol,Crowd);
end