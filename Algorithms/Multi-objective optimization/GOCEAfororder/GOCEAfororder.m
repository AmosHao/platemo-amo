classdef GOCEAfororder< ALGORITHM
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
           pm=1.0/varDim;
           bounds=[Problem.lower;Problem.upper];
           %% Generate random population
           Population = Problem.Initialization();
           pop = Population.decs;
           objVals = Population.objs;
           clustTag=(1:popSize)'; clustName=clustTag; centroid=pop;
           %% Optimization
           gen=1;
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
                           idx=randsample(neigSize,1); parents(1:1,:)=neighborhood(idx,:); 
                       else
                           idx=randsample(globalSize,1); parents(1:1,:)=globalClust(idx,:);
                       end
                   else
                   idx=randsample(globalSize,1); parents(1:1,:)=globalClust(idx,:);
                   end
                   parents(2,:)=currentSol;
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
                   auxPop(i,:)=trialSol; auxVals(i,:)=trialVal;
               end
               auxObjvs=auxVals.objs;
               [auxall,pop,objVals,clustTag,clustName,centroid]=GeOACES_hv_a(auxPop,auxObjvs,pop,objVals,clustTag,clustName,popSize,objDim,varDim,Kmax,auxVals,Population);
               Population=auxall;
               for i = 1:popSize
                   Population(i).add_clustTag = clustTag(i,:);
                   Population(i).add_clustName = clustName(:,:);
               end
           end
        end
    end
end