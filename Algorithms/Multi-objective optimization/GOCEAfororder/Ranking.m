function rk=Ranking(objvs)
% Assign ranks for the population individuals (solutions)
%
%
% Input:
% - objvs - the objective values of the solutions waiting for the rank
% assignment
%
% Output:
% - ranking - the rank assigned for each solution
%
% Author: Zhang Hu, Harbin Institute of Technology
% Last modified: October 27, 2013
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

    %%
    [popSize,objDim]=size(objvs);
    rk=zeros(popSize,1);
    P=(1:popSize)';
    i=1;
    while ~isempty(P)
        [ndfIndex,dfIndex]=NonDominatedFront(objvs(P,1:objDim));
        rk(P(ndfIndex),1)=i;
        P=P(dfIndex);
        i=i+1;
    end
    
end

%%
function [ndfIndex,dfIndex]=NonDominatedFront(objvs)
    %% [ndfIndex, dfIndex] = NonDominatedFront(objvs)
    %
    % Return the indexes of the non-dominated front of the objv_size vectors of
    % M function values contained in pop.
    %
    % IMPORTANT:
    %   Considers Minimization of the objective function values!
    %
    % Input:
    % - objvs               - A matrix of objv_size x M, where M is the number
    %                         of objectives, and objv_size is the number of
    %                         objective function value vectors of the
    %                         solutions.
    %
    % Output:
    % - ndfIndex           - the indexes of the non-dominated front
    % - dfIndex            - the indexes of the solutions that are dominated
    %
    % Author: Zhang Hu, Harbin Institute of Technology
    % Last modified: May 10, 2013
    % Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)
    
    %%
    dfIndex=[];
    ndfIndex=[];
    [popSize,objDim]=size(objvs);
    for i=1:popSize
        v=objvs(i,1:objDim);
        repObjv=v(ones(popSize,1),1:objDim);
        comp1=(objvs<=repObjv);
        comp1=sum(comp1,2);
        if sum(comp1==objDim)>1
            dfIndex=[dfIndex;i];
        else
            ndfIndex=[ndfIndex;i];
        end
    end
    if isempty(ndfIndex)
        ndfIndex=dfIndex;
        dfIndex=[];
    end
end