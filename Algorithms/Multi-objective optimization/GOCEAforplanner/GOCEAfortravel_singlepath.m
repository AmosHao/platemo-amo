classdef GOCEAfortravel_singlepath< ALGORITHM
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
                       if neigSize>n_dots
                           idx=randsample(neigSize,n_dots); 
                           parents(1:n_dots,:)=neighborhood(idx,:); 
                       elseif globalSize>n_dots
                           idx=randsample(globalSize,n_dots); 
                           parents(1:n_dots,:)=globalClust(idx,:);
                       else
                           idx=randsample(popSize,n_dots); parents(1:n_dots,:)=pop(idx,:);
                       end
                   else
                   idx=randsample(popSize,n_dots); parents(1:n_dots,:)=pop(idx,:);
                   end

                   [validPoints]=checknode_singlepath_close(currentSol,varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids);%对当前解包含的航迹点进行碰撞检测
                   if any(~validPoints)
                       rowsToDelete = [rowsToDelete, i];
                   % pop(i,:)=[];Population(i,:)=[];objVals(i,:)=[];clustTag(i,:)=[];
                   zeroIndices = find(~validPoints);
                   
                   % for ind=1:size(zeroIndices,2)
                   %     index=zeroIndices(ind);
                   %     trialSol1=currentSol;
                   %     trialSol2=currentSol;
                   %     if index==1
                   %         trialSol1(:,index)=(currentSol(:,index+1)+currentSol(:,index+2))/2-(currentSol(:,index+1)-currentSol(:,index));
                   %         trialSol1(:,2*index)=(currentSol(:,2*index+1)+currentSol(:,2*index+2))/2-(currentSol(:,2*index+1)-currentSol(:,2*index));
                   %         trialSol1(:,3*index)=(currentSol(:,3*index+1)+currentSol(:,3*index+2))/2-(currentSol(:,3*index+1)-currentSol(:,3*index));
                   %     elseif index==varDim/3
                   %         trialSol1(:,index)=(currentSol(:,index-1)+currentSol(:,index-2))/2-(currentSol(:,index-1)-currentSol(:,index));
                   %         trialSol1(:,2*index)=(currentSol(:,2*index-1)+currentSol(:,2*index-2))/2-(currentSol(:,2*index-1)-currentSol(:,2*index));
                   %         trialSol1(:,3*index)=(currentSol(:,3*index-1)+currentSol(:,3*index-2))/2-(currentSol(:,3*index-1)-currentSol(:,3*index));
                   %     else  
                   %         trialSol1(:,index)=(currentSol(:,index-1)+currentSol(:,index+1))/2;
                   %         trialSol1(:,2*index)=(currentSol(:,2*index-1)+currentSol(:,2*index+1))/2;
                   %         trialSol1(:,3*index)=(currentSol(:,3*index-1)+currentSol(:,3*index+1))/2;
                   %     end
                   %     trialSol2(:,index)=(parents(1,index)+parents(2,index))/2;
                   % end
                   %     trialSol=[trialSol1;trialSol2];
                   for ind=1:size(zeroIndices,2)
                       index=zeroIndices(ind);
                       trialSol1=currentSol;
                       % trialSol2=currentSol;
                       random_number = rand() - 0.5;
                       trialSol1(:,index)=currentSol(:,index)+random_number;
                       trialSol1(:,2*index)=currentSol(:,2*index)+random_number;
                       trialSol1(:,3*index)=currentSol(:,3*index)+random_number;
                       % trialSol2(:,index)=currentSol(:,index)-random_number;
                       % trialSol2(:,2*index)=currentSol(:,2*index)-random_number;
                       % trialSol2(:,3*index)=currentSol(:,3*index)-random_number;
                   end
                   % trialSol2=currentSol;
                   % randomIndices = randperm(n_dots, n_dots);
                   % for g=1:n_dots
                   % trialSol2(:,randomIndices(:,g))=parents(g,randomIndices(:,g));
                   % trialSol2(:,2*randomIndices(:,g))=parents(g,2*randomIndices(:,g));
                   % trialSol2(:,3*randomIndices(:,g))=parents(g,3*randomIndices(:,g));
                   % end

                   %  x=[dot1(1),currentSol(1:n_dots),dot2(1)];
                   %  y=[dot1(2),currentSol(n_dots+1:2*n_dots),dot2(2)];
                   %  z=[dot1(3),currentSol(2*n_dots+1:3*n_dots),dot2(3)];
                   %  for f=1:n_dots
                   %      trialSol_x(:,f)=(x(:,f)+x(:,f+2))/2;
                   %      trialSol_y(:,f)=(y(:,f)+y(:,f+2))/2;
                   %      trialSol_z(:,f)=(z(:,f)+z(:,f+2))/2;
                   %  end
                   % trialSol2=[trialSol_x,trialSol_y,trialSol_z];

                   % trialSol2=currentSol;
                   % random_number = rand() - 0.1;%微小扰动项
                   % trialSol2(:,:)=currentSol(:,:)+random_number;


                   
                   % random_number = 0.99 + (1.01 - 0.99) * rand();%微小扰动项
                   % random_number = rand() - 0.5;
                   % trialSol2(:,:)=currentSol(:,:)+random_number;
                   % trialSol=[trialSol1;trialSol2];
                   trialSol=trialSol1;

                   % else
                   % trialSol=currentSol;
                   % randomIndices = randperm(n_dots, n_dots);
                   % for g=1:n_dots
                   % trialSol(:,randomIndices(:,g))=parents(g,randomIndices(:,g));
                   % trialSol(:,2*randomIndices(:,g))=parents(g,2*randomIndices(:,g));
                   % trialSol(:,3*randomIndices(:,g))=parents(g,3*randomIndices(:,g));
                   % end
                   % % trialSol(:,randomIndices(:,3))=parents(1,randomIndices(:,3));
                   % % trialSol(:,randomIndices(:,4))=parents(2,randomIndices(:,4));
                   % % trialSol(:,randomIndices(:,5))=parents(1,randomIndices(:,5));
                   % % trialSol(:,randomIndices(:,6))=parents(2,randomIndices(:,6));
                   % % count_changedots=ceil(n_dots*0.1);
                   % % randomIndices = randperm(n_dots, count_changedots);
                   % % half_count_changedots=ceil(count_changedots/2);
                   % % for g=1:half_count_changedots
                   % %     trialSol(:,randomIndices(:,g))=parents(1,randomIndices(:,g));
                   % % end
                   % % for h=half_count_changedots+1:count_changedots
                   % %     trialSol(:,randomIndices(:,h))=parents(2,randomIndices(:,h));
                   % % end
                   
                   % trialSol=[trialSol1;trialSol2];
                   else

                   % trialSol=currentSol;
                   % random_number = 0.95 + (1.05 - 0.95) * rand();%微小扰动项
                   % trialSol(:,:)=currentSol(:,:)*random_number;

                   
                   %  % x=[dot1(1),currentSol(1:n_dots),dot2(1)];
                   %  % y=[dot1(2),currentSol(n_dots+1:2*n_dots),dot2(2)];
                   %  % z=[dot1(3),currentSol(2*n_dots+1:3*n_dots),dot2(3)];
                   %  % trialSol_x=[];
                   %  % for f=1:n_dots
                   %  %     trialSol_x(:,f)=(x(:,f)+x(:,f+2))/2;
                   %  %     trialSol_y(:,f)=(y(:,f)+y(:,f+2))/2;
                   %  %     trialSol_z(:,f)=(z(:,f)+z(:,f+2))/2;
                   %  % end
                   %  % trialSol=[trialSol_x,trialSol_y,trialSol_z];
                   % 
                   trialSol=currentSol;
                   randomIndices = randperm(n_dots, n_dots);
                   for g=1:n_dots
                   trialSol(:,randomIndices(:,g))=parents(g,randomIndices(:,g));
                   trialSol(:,2*randomIndices(:,g))=parents(g,2*randomIndices(:,g));
                   trialSol(:,3*randomIndices(:,g))=parents(g,3*randomIndices(:,g));
                   end 
                   end

                   % trialSol2(:,:)=parents(:,:);
                   % random_number = rand() - 0.5;%微小扰动项
                   % trialSol2(:,:)=trialSol2(:,:)+random_number;
                   % trialSol=[trialSol1;trialSol2];
                   % trialVal=Problem.Evaluation(trialSol);
                   % auxPop=[auxPop;trialSol]; auxVals=[auxVals,trialVal];
                   trialVal=Problem.Evaluation(trialSol);
                   auxPop=[auxPop;trialSol]; auxVals=[auxVals,trialVal];  
               end
               pop(rowsToDelete, :)=[];Population(:, rowsToDelete)=[];objVals(rowsToDelete, :)=[];clustTag(rowsToDelete, :)=[];
               % if length(rowsToDelete)>=1
               % trialSol_new=[];
               % a=length(rowsToDelete);
               % popall=[pop;auxPop];
               % select_popall=popall(1:a,:);
               % random_number = rand() - 0.5;%微小扰动项
               % trialSol_new(:,:)=select_popall(:,:)+random_number;
               % trialVal_new=Problem.Evaluation(trialSol_new);
               % auxPop=[auxPop;trialSol_new]; auxVals=[auxVals,trialVal_new]; 
               % end
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