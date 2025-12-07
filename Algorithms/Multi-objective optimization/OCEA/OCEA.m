classdef OCEA < ALGORITHM
% <multi> <real/integer>
% OCEA
% Kmax ---   5 --- The K
% BETA    --- 0.7 --- The BETA
% F ---   0.5 --- The F
% CR    --- 1 --- The CR
    methods
        function main(Algorithm, Problem)
           %% parameter setting
           [Kmax,BETA,F,CR] = Algorithm.ParameterSet(5,0.7,0.5,1);
           popSize = Problem.N;
           % Kmax = 5;
           % BETA = 0.9;
           % F=0.5;CR=1; 
           objDim=Problem.M;
           varDim=Problem.D;
           pm=1.0/varDim;
           bounds=[Problem.lower;Problem.upper];
           %% Generate random population
           Population = Problem.Initialization();
           pop = Population.decs;
           objVals = Population.objs;
           clustTag=(1:popSize)'; counter=[clustTag,ones(popSize,1)]; centroid=pop;
           %% Optimization
           genOCEA=0;
           while Algorithm.NotTerminated(Population)
               genOCEA=genOCEA+1
               for i=1:popSize%遍历整个种群中的个体，为每个个体确定其邻居集合
                   currentSol=pop(i,:); currentTag=clustTag(i);
                   neighborhood=pop(clustTag==currentTag,1:varDim);
                   [~,loc]=ismember(currentSol,neighborhood,'rows');neighborhood(loc,:)=[];                                                   % Determine the neighbouring solutions for the current solution
                   neigSize=size(neighborhood,1);
                   if rand<BETA
                        if neigSize>1
                            idx=randsample(neigSize,2); parents(1:2,1:varDim)=neighborhood(idx,1:varDim);                                      % Select the parents from the neighborhood or the global cluster
                        else
                            idx=randsample(popSize,2); parents(1:2,1:varDim)=pop(idx,1:varDim);
                        end
                   else
                        idx=randsample(popSize,2); parents(1:2,1:varDim)=pop(idx,1:varDim);
                   end
                   parents(3,:)=currentSol;
                   trialSol=DifferentialEvolutionCrossover(parents,bounds,F,CR); trialSol=PolynomialMutation(trialSol,bounds,pm);
                   trialVal=Problem.Evaluation(trialSol);
                   auxVals=trialVal;
                   auxObjvs=trialVal.objs;
                   [auxall,pop,objVals,clustTag,counter,centroid]=OACES_hv_a(trialSol,auxObjvs,pop,objVals,clustTag,counter,centroid,popSize,objDim,varDim,Kmax,auxVals,Population);
               Population=auxall;
               for j = 1:popSize
                   Population(j).add_clustTag = clustTag(j,:);
                   Population(j).add_clustName = counter(:,1);
               end
               end
           end
        end
    end
end