classdef PRMO< ALGORITHM
% <multi> <real/integer>
% GOCEA
% Kmax ---   11 --- The K
% BETA    --- 1 --- The B
% F ---   0.5 --- The F
% CR    --- 1 --- The CR
% maxgen    --- 300 --- The CR
    methods
        function main(Algorithm, Problem)
           %% parameter setting
           [Kmax,BETA,F,CR,maxgen] = Algorithm.ParameterSet(11,1,0.5,1,300);
           popSize = Problem.N; 
           objDim=Problem.M;
           varDim=Problem.D;
           n_dots=varDim/3;
           pm=1.0/varDim;
           bounds=[Problem.lower;Problem.upper];
           % dot1=[500,500,800];dot2=[7000,3000,800];
           dots=Problem.dots;
            dot1=[dots(1,1),dots(1,2),dots(1,3)];
            dot2=[dots(1,4),dots(1,5),dots(1,6)];
            % 读取cubes数据
            %  % 文件名
            % filename = 'city_environment_data_new14_simple.xlsx';
            % 读取长方体数据
            % buildings = xlsread(filename, 'Buildings');
            % numBuildings = size(buildings, 1);
            % 读取圆柱体数据
            % cylinders = xlsread(filename, 'Cylinders');
            % numCylinders = size(cylinders, 1);
            % 读取四棱锥数据
            % pyramids = xlsread(filename, 'Pyramids');
            % numPyramids = size(pyramids, 1);
            % 读取球体数据
            % spheres = xlsread(filename, 'Spheres');
            % numSpheres = size(spheres, 1);
            % 读取禁飞区数据
            % jinfei = xlsread(filename, 'jinfei');
            % numjinfei = size(jinfei, 1);
            % 场景设置
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
           %% Generate random population
           sumb=0;
           while sumb< 2
           Population = Problem.Initialization();
           pop = Population.decs;
           objVals = Population.objs;
           clustTag=(1:popSize)'; clustName=clustTag; centroid=pop;
           a=[];b=[];
                for i = 1:popSize   
                [pop(i,:),validPoints]=checknode_singlepath_close_new(pop(i,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);%对当前解包含的航迹点进行碰撞检测     
                a=[a;validPoints'];
                if any(~validPoints) 
                b=[b;0];
                else   
                b=[b;1];
                end
                end
                sumb=sum(b);
           end
           
           %% Optimization
           gen=0;
           while Algorithm.NotTerminated(Population)
               gen=gen+1
           if gen==1
                a=a;b=b;
           else
               a=auxa;
               b=auxb;
           end
           
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
                   if rnd<BETA*gen/maxgen
                       if neigSize>1
                           idx=randsample(neigSize,2); parents(1:2,:)=neighborhood(idx,:); 
                       else
                           idx=randsample(globalSize,2); parents(1:2,:)=globalClust(idx,:);
                       end
                   else
                   idx=randsample(globalSize,2); parents(1:2,:)=globalClust(idx,:);
                   end

                   
                   % [validPoints]=checknode_singlepath_close(currentSol,varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei);%对当前解包含的航迹点进行碰撞检测
                   validPoints=a(i,:);
                   if any(~validPoints)
                   rowsToDelete = [rowsToDelete, i];
                   zeroIndices = find(~validPoints);
                   for ind=1:size(zeroIndices,2)
                       index=zeroIndices(ind);
                       trialSol1=currentSol;
                       guide_indexs = find(a(:, index) == 1);
                       % guide_indexs_next1 = find(a(:, index+1) == 1);
                       % guide_indexs_next2 = find(a(:, index+1) == 1);
                       if length(guide_indexs)>=1

                       currentdot=[currentSol(:,index),currentSol(:,2*index),currentSol(:,3*index)];
                       dis_guide = zeros(length(guide_indexs), 2);
                       for gg=1:length(guide_indexs)
                       dis_guide(gg,2)=norm([pop(gg,index),pop(gg,2*index),pop(gg,3*index)]-currentdot);
                       dis_guide(gg,1)=guide_indexs(gg,1);
                       end
                       % 找到距离最小的点对应的 guide_index
                        [~, min_dist_idx] = min(dis_guide(:, 2));
                        guide_index = dis_guide(min_dist_idx, 1);

                       % guide_index = guide_indexs(randi(length(guide_indexs)));
                       pop_guide=pop(guide_index,:);
                       trialSol1(:,index)=currentSol(:,index)+(pop_guide(:,index)-currentSol(:,index))*0.1;
                       trialSol1(:,2*index)=currentSol(:,2*index)+(pop_guide(:,2*index)-currentSol(:,2*index))*0.1;
                       trialSol1(:,3*index)=currentSol(:,3*index)+(pop_guide(:,3*index)-currentSol(:,3*index))*0.1;
                       else
                       random_number = rand()*1;
                       trialSol1(:,index)=currentSol(:,index)+random_number;
                       trialSol1(:,2*index)=currentSol(:,2*index)+random_number;
                       trialSol1(:,3*index)=currentSol(:,3*index)+random_number;
                       end
                   end

                   parents(3,:)=currentSol;
                   trialSol2=DifferentialEvolutionCrossover(parents,bounds,F,CR); 
                   trialSol2=PolynomialMutation(trialSol2,bounds,pm);

                   trialSol=[trialSol1;trialSol2];

                   auxPop=[auxPop;trialSol];

                   else
                   parents(3,:)=currentSol;
                   trialSol=DifferentialEvolutionCrossover(parents,bounds,F,CR); 
                   trialSol=PolynomialMutation(trialSol,bounds,pm);
                   auxPop=[auxPop;trialSol]; 

                   end  
                  
               end
               % rowsToDelete3=[rowsToDelete,rowsToDelete2];
               if size(rowsToDelete,2)~=size(pop,1)
               pop(rowsToDelete, :)=[];
               end
               % Population(:, rowsToDelete)=[];objVals(rowsToDelete, :)=[];clustTag(rowsToDelete, :)=[];
               % pop(rowsToDelete3, :)=[];Population(:, rowsToDelete3)=[];objVals(rowsToDelete3, :)=[];clustTag(rowsToDelete3, :)=[];
               % auxObjvs=auxVals.objs;
               selectedSize=Problem.N;

               [auxa,auxb,auxc,auxall,pop,objVals,clustTag,clustName,centroid]=GeOACESfortravel_singlepath_5(auxPop,pop,clustTag,clustName,selectedSize,objDim,varDim,Kmax,Problem,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);
               Population=auxall;
               popSize=size(pop,1);

               for i = 1:popSize
                   Population(i).add_clustTag = clustTag(i,:);
                   Population(i).add_clustName = clustName(:,:);
                   Population(i).validtrait = auxb(i,:);
                   % Population(i).add = auxb(i,:);
                   Population(i).bj = auxc(i,:);

               end

Problem.count_rowstodelete=length(rowsToDelete);

           end
        end
        end
end