classdef simplePRGOCEA< ALGORITHM
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
           %场景设置
        buildings = load('buildings.mat').buildings;
        numBuildings = size(buildings, 1);
        % 读取圆柱体数据
       cylinders = load('cylinders.mat').cylinders;
        numCylinders = size(cylinders, 1);
        % 读取四棱锥数据
        pyramids = load('pyramids.mat').pyramids;
        numPyramids = size(pyramids, 1);
        % 读取球体数据
        spheres = load('spheres.mat').spheres;
        numSpheres = size(spheres, 1);
        % 读取禁飞区数据
        jinfei = load('jinfei.mat').jinfei;
        numjinfei = size(jinfei, 1);
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
    %     % 在 gen=... 时输出数据到 Excel 文件
    %         if gen == 15
    %     writematrix(pop, 'gen15.xlsx', 'Sheet', 1, 'Range', 'A1'); % 输出 popdec
    %     writematrix(clustTag, 'gen15.xlsx', 'Sheet', 2, 'Range', 'A1'); % 输出 clustag
    %     writematrix(clustName, 'gen15.xlsx', 'Sheet', 3, 'Range', 'A1'); % 输出 clustname
    % 
    % end
    % if gen == 100
    %     writematrix(pop, 'gen100.xlsx', 'Sheet', 1, 'Range', 'A1'); % 输出 popdec
    %     writematrix(clustTag, 'gen100.xlsx', 'Sheet', 2, 'Range', 'A1'); % 输出 clustag
    %     writematrix(clustName, 'gen100.xlsx', 'Sheet', 3, 'Range', 'A1'); % 输出 clustname
    % 
    % end
    % 
    %     if gen == 200
    %     writematrix(pop, 'gen200.xlsx', 'Sheet', 1, 'Range', 'A1'); % 输出 popdec
    %     writematrix(clustTag, 'gen200.xlsx', 'Sheet', 2, 'Range', 'A1'); % 输出 clustag
    %     writematrix(clustName, 'gen200.xlsx', 'Sheet', 3, 'Range', 'A1'); % 输出 clustname
    % 
    %     end
    % 
    %         if gen == 300
    %     writematrix(pop, 'gen300.xlsx', 'Sheet', 1, 'Range', 'A1'); % 输出 popdec
    %     writematrix(clustTag, 'gen300.xlsx', 'Sheet', 2, 'Range', 'A1'); % 输出 clustag
    %     writematrix(clustName, 'gen300.xlsx', 'Sheet', 3, 'Range', 'A1'); % 输出 clustname
    % 
    % end

               for i = 1:popSize
                   Population(i).add_clustTag = clustTag(i,:);
                   Population(i).add_clustName = clustName(:,:);
               end
               PopDec=Population.decs;
                varDim=size(PopDec,2);
                for j=1:size(PopDec,1)
                [validPoints]=checknode_singlepath_close(PopDec(j,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei);%对当前解包含的航迹点进行碰撞检测      
               if any(~validPoints) 
               Population(j).validtrait = 0;
               else   
               Population(j).validtrait = 1;
               end
               end
           end
        end
    end
end