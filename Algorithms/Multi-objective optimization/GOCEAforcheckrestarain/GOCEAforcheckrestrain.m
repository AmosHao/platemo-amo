classdef GOCEAforcheckrestrain< ALGORITHM
% <multi> <real/integer>
% GOCEA
% Kmax ---   11 --- The K
% BETA    --- 1 --- The B
% F ---   0.5 --- The F
% CR    --- 1 --- The CR
%maxgen   --- 300 --- The maxgen
    methods
        function main(Algorithm, Problem)
           %% parameter setting
           [Kmax,BETA,F,CR,maxgen] = Algorithm.ParameterSet(11,1,0.5,1,300);
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
                   if rnd<gen/maxgen
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
               [auxall,pop,objVals,clustTag,clustName,centroid]=GeOACES_hv_a(auxPop,auxObjvs,pop,objVals,clustTag,clustName,popSize,objDim,varDim,Kmax,auxVals,Population);
               Population=auxall;


               for i = 1:popSize
    % [validPoints]=checknode_singlepath_close(pop(i,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);%对当前解包含的航迹点进行碰撞检测      
    % a=[a;validPoints'];
    % c=[c;bj'];
    %    if any(~validPoints) 
    %    b=[b;0];
    %    else   
    %    b=[b;1];
    %    end
                   Population(i).add_clustTag = clustTag(i,:);
                   Population(i).add_clustName = clustName(:,:);
                   % Population(i).validtrait = b(i,:);
               end
           end
        end
    end
end