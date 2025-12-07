function [pop,objvs,clustTag,clustName,centroid]=GeOACES_hvBETA(auxPop,auxObjvs,pop,objvs,clustTag,clustName,selectedSize,objDim,varDim,Kmax)
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
pop=[pop;auxPop];tmpSize=size(auxObjvs,1);objvs=[objvs;auxObjvs];%将产生的所有新解加入种群
clustTag=[clustTag;inf(tmpSize,1)];%更新新解类编号为inf
[rk,~]=NDSort(objvs,inf);                                                      % Rank the population                           % Construct the reference point
refPoint=max(objvs,[],1); refPoint=1.2*refPoint; 
auxVals=[]; auxPop=[]; auxTag=[];auxSize=0;                                % Add solutions into the auxiliary population
i=1; FiSize=sum(rk==i);% 按非支配排序编号选取大于N的个体，写入auxVals，auxPop，auxTag
while auxSize<selectedSize
    auxVals=[auxVals;objvs(rk==i,1:objDim)]; auxPop=[auxPop;pop(rk==i,1:varDim)]; auxTag=[auxTag;clustTag(rk==i,1)];auxSize=auxSize+FiSize;
    i=i+1;FiSize=sum(rk==i);
end
FiSize=sum(rk==i-1);

%%hv选择
%
% pop=auxPop;
% objvs=auxVals;
% clustTag=auxTag;
% 
% %%删除具有最小hv的个体，直至pop数等于N
% if size(pop,1)>selectedSize
%     deltaS = CalHV(pop,max(pop,[],1)*1.1,1,10000);
%     while size(pop,1)>selectedSize
%         [~,worst] = min(deltaS);
%         pop(worst,:)=[];
%         objvs(worst,:)=[];
%         clustTag(worst,:)=[];
%         deltaS(worst) = [];
%     end
% end

%%选择具有最大hv的个体，直至pop数等于N
% if size(pop,1)>selectedSize
%     deltaS = CalHV(pop,max(pop,[],1)*1.1,1,10000);
%     %找到 deltaS 中前 selectedSize 个最大值及其索引
%     [maxValues, indices] = maxk(deltaS, selectedSize);
%     %获取对应于 pop 的个体索引
%     selectedIndices = indices;
%     popnew=[];
%     clustTag_new=[];
%     for i=1:length(selectedIndices)
%         index=selectedIndices(i);
%         popnew = [popnew;pop(index,:)];
%         clustTag_new=[clustTag_new;clustTag(index,:)];
%     end
%     pop=popnew;
%     clustTag=clustTag_new;
% end


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
        PSize=PSize-1;
    end
else
    FiObjvs=objvs(rk==i-1,1:objDim); FiPop=pop(rk==i-1,1:varDim); FiTag=clustTag(rk==i-1,1);
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
            FiSize=FiSize-1;
            auxSize=auxSize-1;
        end
    else
        % while auxSize-selectedSize>0
        %     totalHV=HV_test(FiObjvs,refPoint);
        %     %totalHV = CalHV(FiObjvs,max(FiObjvs,[],1)*1.1,1,10000);
        %     fitV=zeros(FiSize,1);
        %     for i=1:FiSize
        %         tmpFiObjvs=FiObjvs;
        %         tmpFiObjvs(i,:)=[];
        %         tmpHV=HV_test(tmpFiObjvs,refPoint);
        %         %tmpHV = CalHV(tmpFiObjvs,max(tmpFiObjvs,[],1)*1.1,1,10000);
        %         fitV(i,:)=totalHV-tmpHV;
        %     end
        %     [~,loc]=min(fitV);
        %     FiObjvs(loc(1,1),:)=[]; FiPop(loc(1,1),:)=[]; FiTag(loc(1,1),:)=[];
        %     FiSize=FiSize-1;
        %     auxSize=auxSize-1;
        % end


        if auxSize-selectedSize>0
            totalHV=HV_test(FiObjvs,refPoint);
            %totalHV = CalHV(FiObjvs,max(FiObjvs,[],1)*1.1,1,10000);
            fitV=zeros(FiSize,1);
            for i=1:FiSize
                tmpFiObjvs=FiObjvs;
                tmpFiObjvs(i,:)=[];
                tmpHV=HV_test(tmpFiObjvs,refPoint);
                %tmpHV = CalHV(tmpFiObjvs,max(tmpFiObjvs,[],1)*1.1,1,10000);
                fitV(i,:)=totalHV-tmpHV;
            end
            [~, sortedIdx] = sort(fitV);
            worst = sortedIdx(1:auxSize-selectedSize);
            %[~,loc]=mink(fitV,auxSize-selectedSize);     
            FiObjvs(worst,:)=[]; FiPop(worst,:)=[]; FiTag(worst,:)=[];
            % FiSize=FiSize-1;
            % auxSize=auxSize-1;
        end


        % %%删除具有最小hv的个体，直至pop数等于N
        % if auxSize-selectedSize>0
        %     deltaS = CalHV(FiPop,max(FiPop,[],1)*1.1,1,10000);
        %     [~, sortedIdx] = sort(deltaS);
        %     worst = sortedIdx(1:auxSize-selectedSize);
        %     % [~,worst] = mink(deltaS,auxSize-selectedSize);
        %     for i = 1:size(worst,2)
        %         FiPop(worst(1,i),:)=[];
        %         FiObjvs(worst(1,i),:)=[];
        %         FiTag(worst(1,i),:)=[];
        %         deltaS(worst(1,i)) = [];
        %     end
        % end




    end
    pop=FiPop; objvs=FiObjvs; clustTag=FiTag;  % The selected population
end




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
pop=pop(clustTag~=inf,1:varDim); objvs=objvs(clustTag~=inf,1:objDim); clustTag=clustTag(clustTag~=inf,1);
for i=1:size(popa,1)
    y=popa(i,1:varDim); yObjv=objvsa(i,1:objDim);
    [pop,objvs,clustTag,clustName,centroid]=OACES_(y,yObjv,pop,objvs,clustTag,clustName,centroid,varDim,Kmax);
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

function [pop,objvs,clustTag,clustName,centroid]=OACES_(y,yObjv,pop,objvs,clustTag,clustName,centroid,varDim,Kmax)
if isempty(clustName)%如果所有先前分好类的个体为空，即全部被删除
    pop=[pop;y]; objvs=[objvs;yObjv]; clustTag=1; clustName=1; centroid=y;
else
    newName=max(clustName)+1;
    pop=[pop;y]; objvs=[objvs;yObjv]; clustTag=[clustTag;max(clustName)+1];
    clustName=[clustName; newName]; centroid=[centroid;y];               % Take the new solution as the centroid of a new cluster
    
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