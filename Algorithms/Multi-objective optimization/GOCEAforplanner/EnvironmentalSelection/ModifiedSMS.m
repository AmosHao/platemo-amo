function [pop,objvs]=ModifiedSMS(pop,objvs,selectedSize,objDim,varDim)
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
    % Rank the population
    rk=paretoRank(objvs);

    % Construct the reference point
    refPoint=max(objvs,[],1);
    refPoint=refPoint+1;

    % Add solutions into the auxiliary population
    auxVals=[];
    auxPop=[];
    auxSize=0;
    i=1;
    FiSize=sum(rk==i);
    while (auxSize+FiSize)<selectedSize
        auxVals=[auxVals;objvs(rk==i,1:objDim)];
        auxPop=[auxPop;pop(rk==i,1:varDim)];
        auxSize=auxSize+FiSize;
        i=i+1;
        FiSize=sum(rk==i);
    end

    if i~=1
        % Solutions in the i-th rank
        FiObjvs=objvs(rk==i,1:objDim);
        FiPop=pop(rk==i,1:varDim);
        pop=[auxPop;FiPop];
        objvs=[auxVals;FiObjvs];
        d=NumberOfDominatingPoints(objvs);
        [~,I]=sort(d,'descend');
        % Remove the solutions with highest d
        objvs(I(1:auxSize+FiSize-selectedSize,1),:)=[];
        pop(I(1:auxSize+FiSize-selectedSize,1),:)=[];
        % The selected population
        
    else
        % Solutions in the i-th rank
        FiObjvs=objvs(rk==i,1:objDim);
        FiPop=pop(rk==i,1:varDim);
        if objDim==2
            for j=1:auxSize+FiSize-selectedSize
                frontObjvs=FiObjvs;
                fitV=zeros(FiSize,1);
                [frontObjvs,IX]=sortrows(frontObjvs,1);
                fitV(IX(1))=(frontObjvs(2,1)-frontObjvs(1,1)).* (refPoint(2)-frontObjvs(1,2));
                fitV(IX(2:FiSize-1))=(frontObjvs(3:FiSize,1)-frontObjvs(2:FiSize-1,1)).* (frontObjvs(1:FiSize-2,2)-frontObjvs(2:FiSize-1,2));
                fitV(IX(FiSize))=(refPoint(1)-frontObjvs(FiSize,1)).*(frontObjvs(FiSize-1,2)-frontObjvs(FiSize,2));
                [~,loc]=min(fitV);
                FiObjvs(loc(1,1),:)=[];
                FiPop(loc(1,1),:)=[];
                FiSize=FiSize-1;
            end
        else
            for j=1:auxSize+FiSize-selectedSize
                totalHV=hv(FiObjvs',refPoint);
                fitV=zeros(FiSize,1);
                for i=1:FiSize
                    tmpObjvs=FiObjvs;
                    tmpObjvs(i,:)=[];
                    tmpHV=hv(tmpObjvs',refPoint);
                    fitV(i,:)=totalHV-tmpHV;
                end
                [~,loc]=min(fitV);
                FiObjvs(loc(1,1),:)=[];
                FiPop(loc(1,1),:)=[];
                FiSize=FiSize-1;
            end
        end
        % The selected population
        pop=[auxPop;FiPop];
        objvs=[auxVals;FiObjvs];
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