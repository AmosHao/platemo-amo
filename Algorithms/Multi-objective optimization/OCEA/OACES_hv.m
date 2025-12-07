function [pop,objvs,clustTag,counter,centroid]=OACES_hv(y,yObjv,pop,objvs,clustTag,counter,centroid,popSize,objDim,varDim,Kmax)
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
objvs=[objvs;yObjv]; pop=[pop;y]; popSize=popSize+1; clustTag=[clustTag;inf];%将产生的所有新解加入种群
[rk,~]=NDSort(objvs,inf);                                                      % Rank the population                           % Construct the reference point
refPoint=max(objvs,[],1); refPoint=1.2*refPoint; 

if max(rk)~=1
        d=NumberOfDominatingPoints(objvs); [~,loc]=max(d);
        if loc~=popSize
        [pop,objvs,clustTag,counter,centroid]=OACES_(y,yObjv,loc,pop,objvs,clustTag,counter,centroid,popSize,objDim,varDim,Kmax); % Delete the solutions with the worst rank and highest d value
        else
         pop(popSize,:)=[]; objvs(popSize,:)=[]; clustTag(popSize,:)=[];
        end
else
    if objDim==2
        frontObjvs=objvs;
        fitV=zeros(popSize,1);
        [frontObjvs,IX]=sortrows(frontObjvs,1);
        fitV(IX(1))=(frontObjvs(2,1)-frontObjvs(1,1)).* (refPoint(2)-frontObjvs(1,2));
        fitV(IX(2:popSize-1))=(frontObjvs(3:popSize,1)-frontObjvs(2:popSize-1,1)).* (frontObjvs(1:popSize-2,2)-frontObjvs(2:popSize-1,2));
        fitV(IX(popSize))=(refPoint(1)-frontObjvs(popSize,1)).*(frontObjvs(popSize-1,2)-frontObjvs(popSize,2));
        [~,loc]=min(fitV);
        if loc(1,1)~=popSize
            [pop,objvs,clustTag,counter,centroid]=OACES_(y,yObjv,loc,pop,objvs,clustTag,counter,centroid,popSize,objDim,varDim,Kmax);
        else
            pop(popSize,:)=[]; objvs(popSize,:)=[]; clustTag(popSize,:)=[];
        end
    else
            totalHV=hv(objvs',refPoint);
            %totalHV = CalHV(FiObjvs,max(FiObjvs,[],1)*1.1,1,10000);
            fitV=zeros(popSize,1);
            for i=1:popSize
                tmpFiObjvs=objvs;
                tmpFiObjvs(i,:)=[];
                tmpHV=hv(tmpFiObjvs',refPoint);
                %tmpHV = CalHV(tmpFiObjvs,max(tmpFiObjvs,[],1)*1.1,1,10000);
                fitV(i,:)=totalHV-tmpHV;
            end
            [~,loc]=min(fitV);
            if loc(1,1)~=popSize
                [pop,objvs,clustTag,counter,centroid]=OACES_(y,yObjv,loc,pop,objvs,clustTag,counter,centroid,popSize,objDim,varDim,Kmax);
            else
                pop(popSize,:)=[]; objvs(popSize,:)=[]; clustTag(popSize,:)=[];
            end
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

function [pop,objvs,clustTag,counter,centroid]=OACES_(y,yObjv,loc,pop,objvs,clustTag,counter,centroid,popSize,objDim,varDim,Kmax)
xTag=clustTag(loc,1);                                                        % Find out that which cluster the loc-th solution locates in
pop(loc,1:varDim)=y; objvs(loc,1:objDim)=yObjv; clustTag(loc,1)=max(counter(:,1))+1; % Remove the worst solution
counter=[counter; [max(counter(:,1))+1,1]]; centroid=[centroid;y];          % Take the new solution as the centroid of a new cluster
pop(popSize,:)=[]; objvs(popSize,:)=[]; clustTag(popSize,:)=[];
% Update the centroid of the cluster that the loc-th solution locates in after removing the loc-th solution
s=sum(clustTag==xTag); 
idx=find(counter(:,1)==xTag);
if s~=0                                                                      % If there still exist the solutions in the cluster x_tag, the centroid and counter will be updated
    counter(idx,2)=counter(idx,2)-1; centroid(idx,1:varDim)=mean(pop(clustTag==xTag,1:varDim),1);
else
    counter(idx,:)=[]; centroid(idx,:)=[];                                   % Otherwise, the corresponding counter and centroid will be deleted
end
nClust=size(counter,1);                                                      % Number of clusters
while nClust>Kmax
    disMatrix=dist(centroid');
    for i=1:nClust
        disMatrix(i,i)=inf;
    end
    minV=min(disMatrix,[],1); [~,ic]=min(minV); [~,ir]=min(disMatrix(:,ic)); % Find two closest clusters
    % Merge the ic-th and ir-th clusters
    newCounter=[max(counter(:,1))+1,counter(ic,2)+counter(ir,2)]; newCentroid=mean([pop(clustTag==counter(ic,1),1:varDim);pop(clustTag==counter(ir,1),1:varDim)],1);
    clustTag(clustTag==counter(ic,1),1)=max(counter(:,1))+1; clustTag(clustTag==counter(ir,1),1)=max(counter(:,1))+1;
    centroid([ic;ir],:)=[]; counter([ic;ir],:)=[];
    counter=[counter;newCounter]; centroid=[centroid;newCentroid];
    nClust=size(counter,1);
end
end
end