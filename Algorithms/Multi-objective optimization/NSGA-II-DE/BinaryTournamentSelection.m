function [selSol,selObjv]=BinaryTournamentSelection(num,pop,objvs)
% Author: Zhang Hu, Harbin Institute of Technology
% Last modified: Nov. 08, 2013
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

    varDim=size(pop,2);objDim=size(objvs,2);
    selSol=rand(num,varDim);selObjv=rand(num,objDim);  % Predefine the selected solutions
    k=0;
    while k<num
        k=k+1;
        [sol,objv,pop,objvs]=Selection(pop,objvs);
        selSol(k,1:varDim)=sol;selObjv(k,1:objDim)=objv;
    end
end

function [selSol,selObjv,pop,objvs]=Selection(pop,objvs)
    [popSize,varDim]=size(pop);
    [~, objDim]=size(objvs);
    if popSize==1
        selSol=pop;selObjv=objvs;
        pop=[];objvs=[];
    else
        rInt=randperm(popSize);% Generate two integers randomly
        randInt1=rInt(1);randInt2=rInt(2);
        % The better solution is selected
        if DominanceComparator(objvs(randInt1,1:objDim),objvs(randInt2,1:objDim))==1
            selSol=pop(randInt1,1:varDim);selObjv=objvs(randInt1,1:objDim);
            pop(randInt1,:)=[];objvs(randInt1,:)=[];
        elseif DominanceComparator(objvs(randInt1,1:objDim),objvs(randInt2,1:objDim))==-1
            selSol=pop(randInt2,1:varDim);selObjv=objvs(randInt2,1:objDim);
            pop(randInt2,:)=[];objvs(randInt2,:)=[];
        elseif rand<0.5
            selSol=pop(randInt1,1:varDim);selObjv=objvs(randInt1,1:objDim);
            pop(randInt1,:)=[];objvs(randInt1,:)=[];
        else
            selSol=pop(randInt2,1:varDim);selObjv=objvs(randInt2,1:objDim);
            pop(randInt2,:)=[];objvs(randInt2,:)=[];
        end
    end
end