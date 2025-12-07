function [pop,objvs]=OriginalSMS(pop,objvs,popSize,nobj)
% Author: Zhang Hu, Harbin Institute of Technology Last modified:
% October 27, 2013 Copyright (C) 2013-2016 by Zhang Hu (e-mail:
% jxzhanghu@126.com)
% This function implements the environmental selection proposed in
% SMS-EMOA. This function is performed for one time, and one solutions is
% removed.

% This method is proposed in SMS-EMOA
refPoint=max(objvs,[],1);
refPoint=refPoint+1;
rk=paretoRank(objvs);
maxRk=max(rk);

if maxRk>1
    d=NumberOfDominatingPoints(objvs);
    % Delete the solutions with the worst rank and highest d value
    [~,loc]=max(d);
    if loc(1,1)~=popSize
        pop(loc(1,1),:)=pop(popSize,:);
        objvs(loc(1,1),:)=objvs(popSize,:);
        pop(popSize,:)=[];
        objvs(popSize,:)=[];
    else
        pop(loc(1),:)=[];
        objvs(loc(1),:)=[];
    end
else
    if nobj==2
        frontObjvs=objvs;
        fitV=zeros(popSize,1);
        [frontObjvs,IX]=sortrows(frontObjvs,1);
        fitV(IX(1))=(frontObjvs(2,1)-frontObjvs(1,1)).* (refPoint(2)-frontObjvs(1,2));
        fitV(IX(2:popSize-1))=(frontObjvs(3:popSize,1)-frontObjvs(2:popSize-1,1)).* (frontObjvs(1:popSize-2,2)-frontObjvs(2:popSize-1,2));
        fitV(IX(popSize))=(refPoint(1)-frontObjvs(popSize,1)).*(frontObjvs(popSize-1,2)-frontObjvs(popSize,2));
        [~,loc]=min(fitV);
        if loc(1,1)~=popSize
            pop(loc(1,1),:)=pop(popSize,:);
            objvs(loc(1,1),:)=objvs(popSize,:);
            pop(popSize,:)=[];
            objvs(popSize,:)=[];
        else
            pop(loc(1),:)=[];
            objvs(loc(1),:)=[];
        end
    else
        totalHV=hv(objvs',refPoint);
        fitV=zeros(popSize,1);
        for i=1:popSize
            tmpObjvs=objvs;
            tmpObjvs(i,:)=[];
            tmpHV=hv(tmpObjvs',refPoint);
            fitV(i,:)=totalHV-tmpHV;
        end
        [~,loc]=min(fitV);
        if loc(1,1)~=popSize
            pop(loc(1,1),:)=pop(popSize,:);
            objvs(loc(1,1),:)=objvs(popSize,:);
            pop(popSize,:)=[];
            objvs(popSize,:)=[];
        else
            pop(loc(1),:)=[];
            objvs(loc(1),:)=[];
        end
    end
end
end

function d=NumberOfDominatingPoints(P)
[PSize,objDim]=size(P);
d=zeros(PSize,1);
% Calculate number of points from P that dominate each solution in S
for i=1:PSize
    sol=P(i,1:objDim);
    repObjv=sol(ones(PSize,1),1:objDim); % copy current individual
    pos=(P<=repObjv);
    pos=sum(pos,2);
    d(i)=sum(pos==objDim)-1;
end
end