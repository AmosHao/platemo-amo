classdef GOCEA_simpleHV < ALGORITHM
% <multi> <real/integer>
% GOCEA
% Kmax ---   9 --- The K
% BETA    --- 0.9 --- The B
% F ---   0.7 --- The F
% CR    --- 1 --- The CR
    methods
        function main(Algorithm, Problem)
           %% parameter setting
           [Kmax,BETA,F,CR] = Algorithm.ParameterSet(9,0.9,0.7,1);
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
           clustTag=(1:popSize)'; clustName=clustTag; centroid=pop;
           gen=0;
           %% Optimization
           while Algorithm.NotTerminated(Population)
               gen=gen+1
               globalClust=[]; %获取全局子种群
               nClust=size(clustName,1);                                              % Number of clusters
               neigSet=cell(nClust,1);
               for i=1:nClust
                   mark=clustTag==clustName(i);%找到当前属于第i类的个体编号
                   neigPop=pop(mark,:);
                   pos=randsample(sum(mark),1);
                   globalClust=[neigPop(pos,:);globalClust];%每个类中随机选一个个体加入到全局子种群
                   if sum(mark)>2
                   neigSet(i,1)={neigPop};
                   end
               end
               globalSize=size(globalClust,1);
               for i=1:popSize%遍历整个种群中的个体，为每个个体确定其邻居集合
                   currentSol=pop(i,:); currentTag=clustTag(i);
                   neighborhood=neigSet{clustName==currentTag,1};
                   if ~isempty(neighborhood)
                   [~,loc]=ismember(currentSol,neighborhood,'rows');neighborhood(loc,:)=[];% Determine the neighbouring solutions for the current solution
                   end %从邻域中删除当前个体
                   neigSize=size(neighborhood,1);
                   % Select the parents from the neighborhood or the global cluster
                   rnd=rand;
                   if rnd<BETA
                       if neigSize>1
                           idx=randsample(neigSize,2); parents(1:2,:)=neighborhood(idx,:); 
                       else
                           idx=randsample(globalSize,2); parents(1:2,:)=globalClust(idx,:);
                       end
                   else
                   idx=randsample(globalSize,2); parents(1:2,:)=globalClust(idx,:);
                   end
                   parents(3,:)=currentSol;
                   trialSol=DifferentialEvolutionCrossover(parents,bounds,F,CR); trialSol=PolynomialMutation(trialSol,bounds,pm);
                   trialVal=Problem.Evaluation(trialSol);
                   auxPop(i,:)=trialSol; auxVals(i,:)=trialVal;
               end
                auxObjvs=auxVals.objs;
               [pop,objVals,clustTag,clustName,centroid]=GeOACES_hv(auxPop,auxObjvs,pop,objVals,clustTag,clustName,popSize,objDim,varDim,Kmax);
               Population=Problem.Evaluation(pop); 
           end
        end
    end
end