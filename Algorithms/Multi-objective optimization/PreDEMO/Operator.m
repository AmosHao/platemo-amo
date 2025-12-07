function Offspring = Operator(Parent1,Parent2,Parent3,FrontNo,Gen,MaxGen,Problem)

%--------------------------------------------------------------------------
% This Operator code consists of Preference-inspired Mutation and Binary
% Crossover and Polynomial Mutation
%% Parameter setting
CR         = 0.15;
F          = 0.5;
proM       = 1;
disM       = 20;
PopObj                = Parent1.objs;
[N1,M]                = size(PopObj);
NF                    = find(FrontNo ==1);
NFs                   = length(NF);

rt                    = sqrt(NFs);
rtm                   = (Gen/MaxGen)^(1/M);
ratio                 = ceil(rt*rtm);


[Pbestind, PbestF]    = FindKneePoints(Parent1(FrontNo==1),ratio);

if isa(Parent1(1),'SOLUTION')
    calObj  = true;
    Parent1 = Parent1.decs;
    Parent2 = Parent2.decs;
    Parent3 = Parent3.decs;
else
    calObj = false;
end

Offspring   = Parent1;
[N,D]   = size(Parent1);

%Problem = PROBLEM.Current();

for i = 1 : N
    Pbestsol     = randi(size(Pbestind,1));
    Pbest        = Pbestind(Pbestsol,:);
    Parent1(i,:) = Parent1(i,:) + F*((Pbest - Parent1(i,:)) + (Parent2(i,:)-Parent3(i,:)));
end

%% Differental evolution
Site = rand(N,D) < CR;
Offspring(Site)  =  Parent1(Site);

%% Polynomial mutation
Lower = repmat(Problem.lower,N,1);
Upper = repmat(Problem.upper,N,1);
Site  = rand(N,D) < proM/D;
mu    = rand(N,D);
temp  = Site & mu<=0.5;
Offspring       = min(max(Offspring,Lower),Upper);
Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
    (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
temp = Site & mu>0.5;
Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
    (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
if calObj
    %Offspring = SOLUTION(Offspring);
    Offspring=Problem.Evaluation(Offspring);
end
end