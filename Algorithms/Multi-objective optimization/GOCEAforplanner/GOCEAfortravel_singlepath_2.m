classdef GOCEAfortravel_singlepath_2< ALGORITHM
% <multi> <real/integer>
% GOCEA
% Kmax ---   15 --- The K
% BETA    --- 0.5 --- The B
% F ---   0.5 --- The F
% CR    --- 1 --- The CR
    methods
        function main(Algorithm, Problem)
           %% parameter setting
           [Kmax,BETA,F,CR] = Algorithm.ParameterSet(15,0.5,0.5,1);
           popSize = Problem.N; 
           objDim=Problem.M;
           varDim=Problem.D;
           n_dots=varDim/3;
           pm=1.0/varDim;
           bounds=[Problem.lower;Problem.upper];
           dot1=[0.2,0.2,0.8];dot2=[2.5,6.8,0.8];
           %% Generate random population
           Population = Problem.Initialization();
           pop = Population.decs;
           objVals = Population.objs;
           clustTag=(1:popSize)'; clustName=clustTag; centroid=pop;
           
           %读取cubes数据
             % 文件名
            filename = 'city_environment_data_6.xlsx';
            % 读取长方体数据
            buildings = xlsread(filename, 'Buildings');
            numBuildings = size(buildings, 1);
            % 读取圆柱体数据
            cylinders = xlsread(filename, 'Cylinders');
            numCylinders = size(cylinders, 1);
            % 读取四棱锥数据
            pyramids = xlsread(filename, 'Pyramids');
            numPyramids = size(pyramids, 1);
            % 读取球体数据
            spheres = xlsread(filename, 'Spheres');
            numSpheres = size(spheres, 1);
           %% Optimization
           while Algorithm.NotTerminated(Population)
               auxPop=[];auxVals=[];
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
               rowsToDelete=[];
               rowsToDelete2=[];
               rowsToDelete3=[];
               for i=1:popSize%遍历整个种群中的个体，为每个个体确定其邻居集合
                   currentSol=pop(i,:); currentTag=clustTag(i);
                   neighborhood=neigSet{clustName==currentTag,1};
                   if ~isempty(neighborhood)
                   [~,loc]=ismember(currentSol,neighborhood,'rows');neighborhood(loc,:)=[];% Determine the neighbouring solutions for the current solution
                   end %从邻域和全局类中删除当前个体
                   % [~,loc]=ismember(currentSol,globalClust,'rows');globalClust(loc,:)=[];
                   % globalSize=size(globalClust,1);
                   neigSize=size(neighborhood,1);
                   % Select the parents from the neighborhood or the global cluster
                   rnd=rand;
                   if rnd<BETA
                       if neigSize>1
                           idx=randsample(neigSize,1); 
                           parents(1,:)=neighborhood(idx,:); 
                       else 
                           [~,loc]=ismember(currentSol,globalClust,'rows');
                           idx=randsample(globalSize,1); 
                           while idx == loc
                                 idx = randsample(globalSize, 1);
                           end
                           parents(1,:)=globalClust(idx,:);
                       end
                   else
                   idx=randsample(globalSize,1);parents(1,:)=globalClust(idx,:);
                   end

                   [validPoints]=checknode_singlepath_close(currentSol,varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids);%对当前解包含的航迹点进行碰撞检测
                   if any(~validPoints)
                   rowsToDelete = [rowsToDelete, i];
                   zeroIndices = find(~validPoints);

                   for ind=1:size(zeroIndices,2)
                       index=zeroIndices(ind);
                       trialSol1=currentSol;
                       random_number = rand() - 0.5;
                       trialSol1(:,index)=currentSol(:,index)+random_number;
                       trialSol1(:,2*index)=currentSol(:,2*index)+random_number;
                       trialSol1(:,3*index)=currentSol(:,3*index)+random_number;

                   end
                   trialSol=trialSol1;
                   trialVal=Problem.Evaluation(trialSol);
                   auxPop=[auxPop;trialSol]; auxVals=[auxVals,trialVal];
                   elseif rnd<0.5
                   % rowsToDelete2 = [rowsToDelete2, i];
                   
                   trialSol=currentSol;
                   mask = logical(mod(1:3*n_dots, 2));
                   trialSol(:, mask) = parents(:, mask);
                   % trialSol1=currentSol;
                   % trialSol2=parents;
                   % randomIndices = randperm(n_dots, 1);
                   % trialSol1(:,randomIndices)=parents(:,randomIndices);
                   % trialSol1(:,2*randomIndices)=parents(:,2*randomIndices);
                   % trialSol1(:,3*randomIndices)=parents(:,3*randomIndices);
                   % trialSol2(:,randomIndices)=currentSol(:,randomIndices);
                   % trialSol2(:,2*randomIndices)=currentSol(:,2*randomIndices);
                   % trialSol2(:,3*randomIndices)=currentSol(:,3*randomIndices);
                   % trialSol=[trialSol1;trialSol2];

                   trialVal=Problem.Evaluation(trialSol);
                   auxPop=[auxPop;trialSol]; auxVals=[auxVals,trialVal];
                   end  
                  
               end
               % rowsToDelete3=[rowsToDelete,rowsToDelete2];
               pop(rowsToDelete, :)=[];Population(:, rowsToDelete)=[];objVals(rowsToDelete, :)=[];clustTag(rowsToDelete, :)=[];
               % pop(rowsToDelete3, :)=[];Population(:, rowsToDelete3)=[];objVals(rowsToDelete3, :)=[];clustTag(rowsToDelete3, :)=[];
               auxObjvs=auxVals.objs;
               selectedSize=Problem.N;
               [auxall,pop,objVals,clustTag,clustName,centroid]=GeOACESfortravel_singlepath(auxPop,auxObjvs,pop,objVals,clustTag,clustName,selectedSize,objDim,varDim,Kmax,auxVals,Population,Problem);
               Population=auxall;
               popSize=size(pop,1);
               for i = 1:popSize
                   Population(i).add_clustTag = clustTag(i,:);
                   Population(i).add_clustName = clustName(:,:);
               end
               % end
Problem.count_rowstodelete=length(rowsToDelete);
           end
        end
        end
end