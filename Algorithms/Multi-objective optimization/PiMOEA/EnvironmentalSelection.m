function [Population,FrontNo,PivotSol] = EnvironmentalSelection(Population,FrontNo,MaxFront,PivotSol,Distance,N)
% This function is used for the environmental selection

% the solutions which are chosen for next generation
Choose = false(1,length(Population));
% all solutions in the first several fronts are chosen
Choose(FrontNo<MaxFront) = 1;
% all the pivot solutions are chosen
Choose(PivotSol) = 1;
if sum(Choose) < N
    % choose some other solutions from the last front if the number of
    % chosen solutions is less than N
    Temp = find(FrontNo==MaxFront & PivotSol==0);
    [null,Rank] = sort(Distance(Temp),'ascend');
    Choose(Temp(Rank(1:(N-sum(Choose))))) = 1;
elseif sum(Choose) > N
    % delete some chosen solutions if the number of chosen solutions is
    % more than N
    Temp = find(FrontNo==MaxFront & PivotSol==1);
    [null,Rank] = sort(Distance(Temp),'descend');
    Choose(Temp(Rank(1:(sum(Choose)-N)))) = 0;
end

Population    = Population(Choose);     % solutions on decision space
FrontNo       = FrontNo(Choose);       % the front number of each solution
PivotSol    = PivotSol(Choose);
end

