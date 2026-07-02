classdef PRMO_fangzhen < ALGORITHM
% <multi> <real/integer>
% 复制自 PRMO.m：仅替换场景为 order_data_fangzhen + buildings_from_osm（详见 planner_fangzhen_scene）
% 运行前请 setenv('PLATEMO_FANGZHEN_ROOT', '你的工程根')，或在 run_17path_PRMO_fangzhen 中已设置
    methods
        function obj = PRMO_fangzhen(varargin)
            % 并行 worker（parfor/spmd）下必须关闭 PlatEMO DefaultOutput（会对 base workspace assignin/save/gui）。
            % outputFcn 只能在构造函数阶段传入（属性对外只读）。
            hasOutFcn = false;
            if ~isempty(varargin)
                isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
                hasOutFcn = any(strcmp(varargin(isStr),'outputFcn'));
            end
            if ~hasOutFcn
                varargin = [varargin, {'outputFcn', @PRMO_fangzhen.silent_output}];
            end
            obj = obj@ALGORITHM(varargin{:});
        end
    end

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
           dots=Problem.dots;
            dot1=[dots(1,1),dots(1,2),dots(1,3)];
            dot2=[dots(1,4),dots(1,5),dots(1,6)];
            repoRoot = getenv('PLATEMO_FANGZHEN_ROOT');
            if isempty(repoRoot)
                repoRoot = 'D:\PlatEMO-master-using';
            end
            addpath(fullfile(repoRoot, 'PlatEMO', 'forOrderNew26', 'Order_Map'));
            scene = planner_fangzhen_scene(repoRoot);
            buildings = scene.buildings;
            % 加速：按起终点的 XY 包围盒（扩 500m）筛选附近建筑，减少碰撞检测数量
            margin = 500;
            minx = min(dot1(1), dot2(1)) - margin;
            maxx = max(dot1(1), dot2(1)) + margin;
            miny = min(dot1(2), dot2(2)) - margin;
            maxy = max(dot1(2), dot2(2)) + margin;
            xcols = 1:3:22;  % 8 个顶点的 x
            ycols = 2:3:23;  % 8 个顶点的 y
            bxmin = min(buildings(:, xcols), [], 2);
            bxmax = max(buildings(:, xcols), [], 2);
            bymin = min(buildings(:, ycols), [], 2);
            bymax = max(buildings(:, ycols), [], 2);
            keep = (bxmax >= minx) & (bxmin <= maxx) & (bymax >= miny) & (bymin <= maxy);
            buildings = buildings(keep, :);
            numBuildings = size(buildings, 1);
            cylinders = scene.cylinders;
            numCylinders = scene.numCylinders;
            pyramids = scene.pyramids;
            numPyramids = scene.numPyramids;
            spheres = scene.spheres;
            numSpheres = scene.numSpheres;
            jinfei = scene.jinfei;
            numjinfei = scene.numjinfei;
           %% Generate random population
           sumb = 0;
           initTry = 0;
           maxInitTries = 30; % 超过该次数仍 sumb<2，则不再循环初始化，直接进入优化
           while sumb < 2
               initTry = initTry + 1;
               Population = Problem.Initialization();
               pop = Population.decs;
               objVals = Population.objs;
               clustTag=(1:popSize)'; clustName=clustTag; centroid=pop;
               a=[];b=[];
               for i = 1:popSize
                   [pop(i,:),validPoints]=checknode_singlepath_close_new(pop(i,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);
                   a=[a;validPoints'];
                   if any(~validPoints)
                       b=[b;0];
                   else
                       b=[b;1];
                   end
               end
               sumb = sum(b);
               if initTry >= maxInitTries
                   break;
               end
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
               globalClust=[]; 
               nClust=size(clustName,1);
               neigSet=cell(nClust,1);
               for i=1:nClust
                   mark=clustTag==clustName(i);
                   neigPop=pop(mark,:);
                   pos=randsample(sum(mark),1);
                   globalClust=[neigPop(pos,:);globalClust];
                   if sum(mark)>2
                   neigSet(i,1)={neigPop};
                   end
               end
               globalSize=size(globalClust,1);
               rowsToDelete=[];
               rowsToDelete2=[];
               rowsToDelete3=[];

               for i=1:popSize
                   currentSol=pop(i,:); currentTag=clustTag(i);
                   neighborhood=neigSet{clustName==currentTag,1};
                   if ~isempty(neighborhood)
                   [tf, loc] = ismember(currentSol, neighborhood, 'rows');
                   if tf && loc > 0
                       neighborhood(loc, :) = [];
                   else
                       % 浮点误差会导致 ismember 失败(loc=0)，退化为删掉与 currentSol 最近的行
                       d2 = sum((neighborhood - currentSol).^2, 2);
                       [~, imin] = min(d2);
                       neighborhood(imin, :) = [];
                   end
                   end
                   neigSize=size(neighborhood,1);
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
                   validPoints=a(i,:);
                   if any(~validPoints)
                   rowsToDelete = [rowsToDelete, i];
                   zeroIndices = find(~validPoints);
                   for ind=1:size(zeroIndices,2)
                       index=zeroIndices(ind);
                       trialSol1=currentSol;
                       guide_indexs = find(a(:, index) == 1);
                       if length(guide_indexs)>=1
                       currentdot=[currentSol(:,index),currentSol(:,2*index),currentSol(:,3*index)];
                       dis_guide = zeros(length(guide_indexs), 2);
                       for gg=1:length(guide_indexs)
                       dis_guide(gg,2)=norm([pop(gg,index),pop(gg,2*index),pop(gg,3*index)]-currentdot);
                       dis_guide(gg,1)=guide_indexs(gg,1);
                       end
                        [~, min_dist_idx] = min(dis_guide(:, 2));
                        guide_index = dis_guide(min_dist_idx, 1);
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
               if size(rowsToDelete,2)~=size(pop,1)
               pop(rowsToDelete, :)=[];
               end
               selectedSize=Problem.N;
               [auxa,auxb,auxc,auxall,pop,objVals,clustTag,clustName,centroid]=GeOACESfortravel_singlepath_5(auxPop,pop,clustTag,clustName,selectedSize,objDim,varDim,Kmax,Problem,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);
               Population=auxall;
               popSize=size(pop,1);
               for i = 1:popSize
                   Population(i).add_clustTag = clustTag(i,:);
                   Population(i).add_clustName = clustName(:,:);
                   Population(i).validtrait = auxb(i,:);
                   Population(i).bj = auxc(i,:);
               end
Problem.count_rowstodelete=length(rowsToDelete);
           end
        end
    end
    methods(Static, Access = private)
        function silent_output(Algorithm, Problem)
        % 静默输出：仅打印一行进度，不做 GUI / 不写 Data/*.mat / 不 assignin base
            fprintf('%s on %d-objective %d-variable %s (%6.2f%%), %.2fs passed...\n', ...
                class(Algorithm), Problem.M, Problem.D, class(Problem), ...
                Problem.FE/Problem.maxFE*100, Algorithm.metric.runtime);
        end
    end
end
