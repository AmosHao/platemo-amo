function [Population,FrontNo] = EnvironmentalSelection(Population,Offspring,N)
% The environmental selection of Pre-DEMO
% Select the offsprings dominating their corresponding parents firstly,
% and then select the offspring population non-dominated with the
%% Select by constraint-domination
PopObj    = Population.objs;
PopCon    = Population.cons;
feasibleP = all(PopCon<=0,2);
OffObj    = Offspring.objs;
OffCon    = Offspring.cons;
feasibleO = all(OffCon<=0,2);
% The offsprings which can replace its parent
updated = ~feasibleP&feasibleO  | ...
    ~feasibleP&~feasibleO & all(PopCon>=OffCon,2) | ...
    feasibleP&feasibleO   & all(PopObj>=OffObj,2);
% The offsprings which can add to the population
selected = feasibleP&feasibleO & any(PopObj<OffObj,2) & any(PopObj>OffObj,2);
% Update the population
Population(updated) = Offspring(updated);
Population          = [Population,Offspring(selected)];

%% Non-dominated sorting wrt the actual objectives
%% Selection among feasible solutions
PopObj            = Population.objs;
[FrontNo,MaxFNo]  = NDSort(PopObj,N);
%% Calculate the crowding distance of each solution
Distance   = CalFitness(PopObj,FrontNo, MaxFNo); 
Fitness    = [FrontNo' Distance'];
[~,Rank]   = sortrows(Fitness);
NSel       = Rank(1:N);
Next(NSel) = true;
%% Population for next generation
Population = Population(Next);
FrontNo    = FrontNo(Next); 


end