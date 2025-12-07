function [auxa,auxb,auxc,auxall,pop,objvs,clustTag,clustName,centroid]=GeOACESfortravel_singlepath_5(auxPop,pop,clustTag,clustName,selectedSize,objDim,varDim,Kmax,Problem,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2,gen)
% Author: Zhang Hu, Harbin Institute of Technology Last modified:
% October 27, 2013 Copyright (C) 2013-2016 by Zhang Hu (e-mail:
% jxzhanghu@126.com)
% The function implements selecting selectedSize solutions from the
% population based on the modified S-metric environmental selection approach
% proposed in SMS-EMOA. Different from the original approach, this approach
% is performed after that the whole population solutions have produced
% their offspring solutions, and all the offspring solutions are combined
% with their parents to form a big population, the approach is conducted on
% the big population to choose selectedSize better solutions.
% INPUT:
% pop: population to be selected
% objvs: objective values of the population
% selectedSize: predefined size of the solutions required to be selected
% objDim: objective dimension
% varDim: variable dimension
% OUTPUT:
% pop: selected population
% objvs: objective values of selected population

% This method is proposed in SMS-EMOA
popinitial=[pop;auxPop];tmpSize=size(auxPop,1);
% objvs=[objvs;auxObjvs];%将产生的所有新解加入种群
% auxVals=auxVals';
% allVals=[Population,auxVals];
clustTag=[clustTag;inf(tmpSize,1)];%更新新解类编号为inf

% [validPoints,pop]=checknode(pop);
% allVals=Problem.Evaluation(pop);
% objvs=allVals.objs; 
a=[];b=[];c=[];
for i = 1:size(popinitial,1)           
    [popinitial(i,:),validPoints,bj]=checknode_singlepath_close_newbj(popinitial(i,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);%对当前解包含的航迹点进行碰撞检测
              % Population(i).dec = pop(i,:);       
    a=[a;validPoints'];
    c=[c;bj'];
       if any(~validPoints) 
       b=[b;0];
            if gen>=0
                zeroIndices = find(~validPoints);
                n=varDim/4;
                popinitial(i,zeroIndices+3*n)=popinitial(i,zeroIndices+3*n)-0.2;
            end
       else   
       b=[b;1];
       end

end
allVals = Problem.Evaluation(popinitial);
objvs = allVals.objs;

% objcheck=objvs;
% deccheck=allVals.decs;
% % validtraitcheck=[];
% % for kk=1:size(deccheck,1)
% % validtraitcheck(kk,:)=allVals(kk).validtrait;
% % end
% validtraitcheck=b;
% check2(objcheck,deccheck,validtraitcheck)

popvalid=popinitial(b==1,:);
popunvalid=popinitial(b==0,:);
objvsvalid=objvs(b==1,:);
objvsunvalid=objvs(b==0,:);
allValsvalid=allVals(:,b==1);
allValsunvalid=allVals(:,b==0);
clustTagvalid=clustTag(b==1,:);
clustTagunvalid=clustTag(b==0,:);
avalid=a(b==1,:);
aunvalid=a(b==0,:);
bvalid=b(b==1,:);
bunvalid=b(b==0,:);
cvalid=c(b==1,:);
cunvalid=c(b==0,:);

[rk,~]=NDSort(objvs,inf);                                                      % Rank the population                           % Construct the reference point
refPoint=max(objvs,[],1); refPoint=1.2*refPoint; 
auxVals=[]; auxPop=[]; auxTag=[];auxSize=0; auxall=[];  auxa=[];auxb=[];   auxc=[];                          % Add solutions into the auxiliary population
i=1; FiSize=sum(rk==i);% 按非支配排序编号选取大于N的个体，写入auxVals，auxPop，auxTag
while auxSize<selectedSize
    auxVals=[auxVals;objvs(rk==i,1:objDim)];
    auxPop=[auxPop;popinitial(rk==i,1:varDim)];
    auxTag=[auxTag;clustTag(rk==i,1)];
    auxSize=auxSize+FiSize;
    auxall=[auxall,allVals(rk==i)];
    auxa=[auxa;a(rk==i,:)];
    auxb=[auxb;b(rk==i,:)];
    auxc=[auxc;c(rk==i,:)];
    i=i+1;FiSize=sum(rk==i);
end
FiSize=sum(rk==i-1);


%原始方法
if i-1~=1
    % FiObjvs=objvs(rk==i,1:objDim); FiPop=pop(rk==i,1:varDim); FiTag=clustTag(rk==i,1); % Solutions in the i-th rank
    % pop=[auxPop;FiPop]; objvs=[auxVals;FiObjvs];clustTag=[auxTag;FiTag];
    %pSize=auxSize+FiSize;
    pop=auxPop;
    objvs=auxVals;
    clustTag=auxTag;
    PSize = auxSize;
    while PSize-selectedSize>0
        d=NumberOfDominatingPoints(objvs); [~,loc]=max(d);
        objvs(loc(1,1),:)=[]; pop(loc(1,1),:)=[]; clustTag(loc(1,1),:)=[];% Remove the solutions with highest d
        auxall(loc(1,1))=[];
            auxa(loc(1,1),:)=[];
            auxb(loc(1,1),:)=[];
            auxc(loc(1,1),:)=[];
        PSize=PSize-1;
    end
else
    % FiObjvs=objvs(rk==i-1,1:objDim); FiPop=pop(rk==i-1,1:varDim); FiTag=clustTag(rk==i-1,1);
    FiObjvs=auxVals; FiPop=auxPop; FiTag=auxTag;
    if objDim==2
        while auxSize-selectedSize>0
            frontObjvs=FiObjvs;
            fitV=zeros(FiSize,1);
            [frontObjvs,IX]=sortrows(frontObjvs,1);
            fitV(IX(1))=(frontObjvs(2,1)-frontObjvs(1,1)).* (refPoint(1,2)-frontObjvs(1,2));
            fitV(IX(2:FiSize-1))=(frontObjvs(3:FiSize,1)-frontObjvs(2:FiSize-1,1)).* (frontObjvs(1:FiSize-2,2)-frontObjvs(2:FiSize-1,2));
            fitV(IX(FiSize))=(refPoint(1,1)-frontObjvs(FiSize,1)).*(frontObjvs(FiSize-1,2)-frontObjvs(FiSize,2));
            [~,loc]=min(fitV);
            FiObjvs(loc(1,1),:)=[]; FiPop(loc(1,1),:)=[]; FiTag(loc(1,1),:)=[];
            auxall(loc(1,1))=[];
            auxa(loc(1,1),:)=[];
            auxb(loc(1,1),:)=[];
            auxc(loc(1,1),:)=[];
            FiSize=FiSize-1;
            auxSize=auxSize-1;
        end
    else
tic;
        while auxSize-selectedSize>0
            totalHV=hvtest(FiObjvs',refPoint);
            %totalHV = CalHV(FiObjvs,max(FiObjvs,[],1)*1.1,1,10000);
            fitV=zeros(FiSize,1);
            for i=1:FiSize
                tmpFiObjvs=FiObjvs;
                tmpFiObjvs(i,:)=[];
                tmpHV=hvtest(tmpFiObjvs',refPoint);
                %tmpHV = CalHV(tmpFiObjvs,max(tmpFiObjvs,[],1)*1.1,1,10000);
                fitV(i,:)=totalHV-tmpHV;
            end     
            [~,loc]=min(fitV);
            FiObjvs(loc(1,1),:)=[]; FiPop(loc(1,1),:)=[]; FiTag(loc(1,1),:)=[];
            auxall(loc(1,1))=[];
            auxa(loc(1,1),:)=[];
            auxb(loc(1,1),:)=[];
            auxc(loc(1,1),:)=[];
            FiSize=FiSize-1;
            auxSize=auxSize-1;
        end
    end
    toc;
    pop=FiPop; objvs=FiObjvs; clustTag=FiTag;  % The selected population
end

if sum(auxb)<selectedSize*0.3
   [rkvalid,~]=NDSort(objvsvalid,inf); 
   [rkunvalid,~]=NDSort(objvsunvalid,inf);
   auxVals=[]; auxPop=[]; auxTag=[];auxSize=0; auxall=[];auxa=[];auxb=[];auxc=[];                             % Add solutions into the auxiliary population
i=1; FiSizevalid=sum(rkvalid==i);FiSizeunvalid=sum(rkunvalid==i);
FiSize=FiSizevalid+FiSizeunvalid;% 按非支配排序编号选取大于N的个体，写入auxVals，auxPop，auxTag
while auxSize<selectedSize
    auxVals=[auxVals;objvsvalid(rkvalid==i,1:objDim)];
    auxVals=[auxVals;objvsunvalid(rkunvalid==i,1:objDim)];
    auxPop=[auxPop;popvalid(rkvalid==i,1:varDim)];
    auxPop=[auxPop;popunvalid(rkunvalid==i,1:varDim)];
    auxTag=[auxTag;clustTagvalid(rkvalid==i,1)];
    auxTag=[auxTag;clustTagunvalid(rkunvalid==i,1)];
    auxSize=auxSize+FiSize;
    auxall=[auxall,allValsvalid(rkvalid==i)];
    auxall=[auxall,allValsunvalid(rkunvalid==i)];
    auxa=[auxa;avalid(rkvalid==i,:)];
    auxa=[auxa;aunvalid(rkunvalid==i,:)];
    auxb=[auxb;bvalid(rkvalid==i,:)];
    auxb=[auxb;bunvalid(rkunvalid==i,:)];
    auxc=[auxc;cvalid(rkvalid==i,:)];
    auxc=[auxc;cunvalid(rkunvalid==i,:)];
    i=i+1;FiSizevalid=sum(rkvalid==i);FiSizeunvalid=sum(rkunvalid==i);
    FiSize=FiSizevalid+FiSizeunvalid;
end
FiSizevalid=sum(rkvalid==i-1);FiSizeunvalid=sum(rkunvalid==i-1);
    FiSize=FiSizevalid+FiSizeunvalid;
%原始方法
if i-1~=1
    % FiObjvs=objvs(rk==i,1:objDim); FiPop=pop(rk==i,1:varDim); FiTag=clustTag(rk==i,1); % Solutions in the i-th rank
    % pop=[auxPop;FiPop]; objvs=[auxVals;FiObjvs];clustTag=[auxTag;FiTag];
    %pSize=auxSize+FiSize;
    pop=auxPop;
    objvs=auxVals;
    clustTag=auxTag;
    PSize = auxSize;
    while PSize-selectedSize>0
        d=NumberOfDominatingPoints(objvs); 
        % [~,loc]=max(d);
        % 找到 auxb 中等于 0 的位置
        zero_indices = find(auxb == 0);
        % 找到这些位置中 d 取最大值的索引
        [~, max_idx] = max(d(zero_indices));
        % 计算最终的索引位置
        loc = zero_indices(max_idx);
        objvs(loc(1,1),:)=[]; pop(loc(1,1),:)=[]; clustTag(loc(1,1),:)=[];% Remove the solutions with highest d
        auxall(loc(1,1))=[];
            auxa(loc(1,1),:)=[];
            auxb(loc(1,1),:)=[];
            auxc(loc(1,1),:)=[];
        PSize=PSize-1;
    end
else
    % FiObjvs=objvs(rk==i-1,1:objDim); FiPop=pop(rk==i-1,1:varDim); FiTag=clustTag(rk==i-1,1);
    FiObjvs=auxVals; FiPop=auxPop; FiTag=auxTag;
    if objDim==2
        while auxSize-selectedSize>0
            frontObjvs=FiObjvs;
            fitV=zeros(FiSize,1);
            [frontObjvs,IX]=sortrows(frontObjvs,1);
            fitV(IX(1))=(frontObjvs(2,1)-frontObjvs(1,1)).* (refPoint(1,2)-frontObjvs(1,2));
            fitV(IX(2:FiSize-1))=(frontObjvs(3:FiSize,1)-frontObjvs(2:FiSize-1,1)).* (frontObjvs(1:FiSize-2,2)-frontObjvs(2:FiSize-1,2));
            fitV(IX(FiSize))=(refPoint(1,1)-frontObjvs(FiSize,1)).*(frontObjvs(FiSize-1,2)-frontObjvs(FiSize,2));
            % [~,loc]=min(fitV);
            % 找到 auxb 中等于 0 的位置
            zero_indices = find(auxb == 0);
            % 找到这些位置中 d 取最大值的索引
            [~, min_idx] = min(fitV(zero_indices));
            % 计算最终的索引位置
            loc = zero_indices(min_idx);
            FiObjvs(loc(1,1),:)=[]; FiPop(loc(1,1),:)=[]; FiTag(loc(1,1),:)=[];
            auxall(loc(1,1))=[];
            auxa(loc(1,1),:)=[];
            auxb(loc(1,1),:)=[];
            auxc(loc(1,1),:)=[];
            FiSize=FiSize-1;
            auxSize=auxSize-1;
        end
    else
tic;
        while auxSize-selectedSize>0
            totalHV=hvtest(FiObjvs',refPoint);
            %totalHV = CalHV(FiObjvs,max(FiObjvs,[],1)*1.1,1,10000);
            fitV=zeros(FiSize,1);
            for i=1:FiSize
                tmpFiObjvs=FiObjvs;
                tmpFiObjvs(i,:)=[];
                tmpHV=hvtest(tmpFiObjvs',refPoint);
                %tmpHV = CalHV(tmpFiObjvs,max(tmpFiObjvs,[],1)*1.1,1,10000);
                fitV(i,:)=totalHV-tmpHV;
            end     
            % [~,loc]=min(fitV);
            % 找到 auxb 中等于 0 的位置
            zero_indices = find(auxb == 0);
            % 找到这些位置中 d 取最大值的索引
            [~, min_idx] = min(fitV(zero_indices));
            % 计算最终的索引位置
            loc = zero_indices(min_idx);
            FiObjvs(loc(1,1),:)=[]; FiPop(loc(1,1),:)=[]; FiTag(loc(1,1),:)=[];
            auxall(loc(1,1))=[];
            auxa(loc(1,1),:)=[];
            auxb(loc(1,1),:)=[];
            auxc(loc(1,1),:)=[];
            FiSize=FiSize-1;
            auxSize=auxSize-1;
        end
    end
    toc;
    pop=FiPop; objvs=FiObjvs; clustTag=FiTag;  % The selected population
end

end

% [validPoints,pop]=checknode(pop);
% auxall=Problem.Evaluation(pop);
% objvs=auxall.objs;


K=size(clustName,1);%对已经分类的个体重新写入类编号和计算类中心
centroidPrime=[]; clustNamePrime=[];
for i=1:K
    if ismember(clustName(i),clustTag)
        clustNamePrime=[clustNamePrime;clustName(i,1)];
        centroidPrime=[centroidPrime;mean(pop(clustTag==clustName(i),1:varDim),1)];
    end
end
clustName=clustNamePrime; centroid=centroidPrime;
popa=pop(clustTag==inf,1:varDim); objvsa=objvs(clustTag==inf,1:objDim);
auxallpopa=auxall(:,clustTag==inf);
auxapopa=auxa(clustTag==inf,:);
auxbpopa=auxb(clustTag==inf,:);
auxcpopa=auxc(clustTag==inf,:);
pop=pop(clustTag~=inf,1:varDim); 
objvs=objvs(clustTag~=inf,1:objDim);
auxall=auxall(:,clustTag~=inf); 
auxa=auxa(clustTag~=inf,:);
auxb=auxb(clustTag~=inf,:);
auxc=auxc(clustTag~=inf,:);
clustTag=clustTag(clustTag~=inf,1);%%%%%%%%%%%错误修改标记
for i=1:size(popa,1)
    y=popa(i,1:varDim); yObjv=objvsa(i,1:objDim);ypopulation=auxallpopa(:,i);
    ya=auxapopa(i,:);
    yb=auxbpopa(i,:);
    yc=auxcpopa(i,:);
    [auxa,auxb,auxc,auxall,pop,objvs,clustTag,clustName,centroid]=OACES_(auxall,y,yObjv,ypopulation,ya,yb,yc,pop,objvs,clustTag,auxa,auxb,auxc,clustName,centroid,varDim,Kmax);
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

function [auxa,auxb,auxc,auxall,pop,objvs,clustTag,clustName,centroid]=OACES_(auxall,y,yObjv,ypopulation,ya,yb,yc,pop,objvs,clustTag,auxa,auxb,auxc,clustName,centroid,varDim,Kmax)
if isempty(clustName)%如果所有先前分好类的个体为空，即全部被删除
    pop=[pop;y]; objvs=[objvs;yObjv];auxall=[auxall,ypopulation]; 
    auxa=[auxa;ya];auxb=[auxb;yb];auxc=[auxc;yc];
    clustTag=1; clustName=1; centroid=y;
else
    newName=max(clustName)+1;
    pop=[pop;y]; objvs=[objvs;yObjv]; clustTag=[clustTag;max(clustName)+1];
    clustName=[clustName; newName]; centroid=[centroid;y];               % Take the new solution as the centroid of a new cluster
    auxall=[auxall,ypopulation];
    auxa=[auxa;ya];auxb=[auxb;yb];auxc=[auxc;yc];
    nClust=size(clustName,1);                                              % Number of clusters
    while nClust>Kmax
        distances=pdist(centroid);
        disMatrix = squareform(distances);
        for i=1:nClust
            disMatrix(i,i)=inf;
        end
        minV=min(disMatrix,[],1); [~,ic]=min(minV); [~,ir]=min(disMatrix(:,ic)); % Find two closest clusters
        % Merge the ic-th and ir-th clusters
        newName=newName+1; newCentroid=mean([pop(clustTag==clustName(ic,1),1:varDim);pop(clustTag==clustName(ir,1),1:varDim)],1);
        clustTag(clustTag==clustName(ic,1),1)=newName; clustTag(clustTag==clustName(ir,1),1)=newName;
        centroid([ic;ir],:)=[]; clustName([ic;ir],:)=[];
        clustName=[clustName;newName]; centroid=[centroid;newCentroid];
        nClust=size(clustName,1);
    end
end

end