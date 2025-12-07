function rankCD=RankCrowding(objvs)
    % Apply non-dominated sorting on the objvs. Assign a rank for each
    % solution and calculate the crowding distance for the solutions in the
    % same ranks.
    %
    % Author: Zhang Hu, Harbin Institute of Technology
    % Last modified: October 27, 2013
    % Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

    %%
    [popSize,objDim]=size(objvs);
    rankCD=zeros(popSize,2);
    P=(1:popSize)';
    i=1;
    while ~isempty(P)
        [ndfIndex,dfIndex]=NonDominatedFront(objvs(P,1:objDim));
        rankCD(P(ndfIndex),1)=i;
        rankCD(P(ndfIndex),2)=Crowding(objvs(P(ndfIndex),1:objDim));
        P=P(dfIndex);
        i=i+1;
    end
end

%%
function [ndfIndex,dfIndex]=NonDominatedFront(objvs)
% Author: Zhang Hu, Harbin Institute of Technology
% Last modified: May 10, 2013
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)
    
    dfIndex=[];ndfIndex=[];
    [popSize,objDim]=size(objvs);
    for i=1:popSize
        v=objvs(i,1:objDim);
        repObjv=v(ones(popSize,1),1:objDim);
        comResults1=(objvs<=repObjv);
        comResults1=sum(comResults1,2);
        if sum(comResults1==objDim)>1
            dfIndex=[dfIndex;i];
        else
            ndfIndex=[ndfIndex;i];
        end
    end
    if isempty(ndfIndex)
        ndfIndex=dfIndex;dfIndex=[];
    end
end

function crowdDis=Crowding(objvs)
% Author: Author: Zhang Hu, Harbin Institute of Technology
% Last modified: October 27, 2013
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

    %%
    [popSize,objDim]=size(objvs);
    if popSize==0
        error('The input is empty!');
    elseif popSize==1
        crowdDis=inf;
    elseif popSize==2
        crowdDis=[inf;inf];
    else
        crowdDis=zeros(popSize,1);
        for j=1:objDim
            if max(objvs(1:popSize,j))==min(objvs(1:popSize,j))
                continue;
            else
                [~,sortIdx]=sort(objvs(1:popSize,j));
                crowdDis(sortIdx(1),1)=inf;
                crowdDis(sortIdx(end),1)=inf;
                for k=2:popSize-1
                    crowdDis(sortIdx(k),1)=crowdDis(sortIdx(k),1)+((objvs(sortIdx(k+1),j)-objvs(sortIdx(k-1),j))/(objvs(sortIdx(end),j)-objvs(sortIdx(1),j)));
                end
            end
        end
    end
end