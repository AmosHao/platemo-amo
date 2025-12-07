classdef GOCEAfororder_1< ALGORITHM
% <multi> <real/integer>
% GOCEA
% Kmax ---   11 --- The K
% BETA    --- 0.8 --- The B
% F ---   0.5 --- The F
% CR    --- 1 --- The CR
    methods
        function main(Algorithm, Problem)
           %% parameter setting
           [Kmax,BETA,F,CR] = Algorithm.ParameterSet(11,0.8,0.5,1);
           popSize = Problem.N; 
           objDim=Problem.M;
           varDim=Problem.D;
           % pm=1.0/varDim;
           % bounds=[Problem.lower;Problem.upper];
           %% Generate random population
           Population = Problem.Initialization();
           pop = Population.decs;
           objvs = Population.objs;
           % clustTag=(1:popSize)'; clustName=clustTag; centroid=pop;
           %% Optimization
           gen=1;
           while Algorithm.NotTerminated(Population)
               gen=gen+1
               idx=randsample(popSize,2); 
               parents=pop(idx,:);
               trialSol=pop(1,:);
               while any(all(trialSol == pop, 2))
                   % trialSol=DifferentialEvolutionCrossover(parents,bounds,F,CR); 
                   [child1, child2]=genetic_crossover(parents(1,:),parents(2,:));
                   % trialSol=PolynomialMutation(trialSol,bounds,pm);
                   childs=[child1;child2];
                   index=randi([1,2]);
                   child=childs(index,:);
                   [child]=bianyi(child);
                   trialSol=child;
                end
                   trialVal=Problem.Evaluation(trialSol);
                   auxPop=trialSol; 
                   auxVals=trialVal;
                   auxObjvs=auxVals.objs;
                   pop=[pop;auxPop];
                   objvs=[objvs;auxObjvs];%将产生的所有新解加入种群
                   auxVals=auxVals';
                   allVals=[Population,auxVals];
                   [rk,~]=NDSort(objvs,inf);
                   refPoint=max(objvs,[],1); refPoint=1.2*refPoint;
if max(rk)~=1
        d=NumberOfDominatingPoints(objvs); [~,loc]=max(d);
        objvs(loc(1,1),:)=[]; pop(loc(1,1),:)=[]; 
        allVals(loc(1,1))=[];
        Population=allVals;
else
    FiSize=sum(rk==max(rk));
    FiObjvs=objvs(rk==max(rk),1:objDim); FiPop=pop(rk==max(rk),1:varDim); 
            frontObjvs=FiObjvs;
            fitV=zeros(FiSize,1);
            [frontObjvs,IX]=sortrows(frontObjvs,1);
            fitV(IX(1))=(frontObjvs(2,1)-frontObjvs(1,1)).* (refPoint(1,2)-frontObjvs(1,2));
            fitV(IX(2:FiSize-1))=(frontObjvs(3:FiSize,1)-frontObjvs(2:FiSize-1,1)).* (frontObjvs(1:FiSize-2,2)-frontObjvs(2:FiSize-1,2));
            fitV(IX(FiSize))=(refPoint(1,1)-frontObjvs(FiSize,1)).*(frontObjvs(FiSize-1,2)-frontObjvs(FiSize,2));
            [~,loc]=min(fitV);
            FiObjvs(loc(1,1),:)=[]; FiPop(loc(1,1),:)=[];
            allVals(loc(1,1))=[];
            objvs=FiObjvs;
            pop=FiPop;
            Population=allVals;
end
           end
        end
    end
end
function d=NumberOfDominatingPoints(P)
[PSize,objDim]=size(P);
d=zeros(PSize,1);
% Calculate number of points from P that dominate each solution in S
for i=1:PSize
    sol=P(i,1:objDim);
    repObjv=sol(ones(PSize,1),1:objDim);                                   % copy current individual
    pos=(P<=repObjv);
    pos=sum(pos,2);
    d(i)=sum(pos==objDim)-1;
end
end