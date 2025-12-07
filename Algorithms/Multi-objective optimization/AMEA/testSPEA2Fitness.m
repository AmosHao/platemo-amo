function rawFit=testSPEA2Fitness(objVals)
% This function implements the fitness assignment depicted in SPEA2.
% Different from our usual concept, in this program, the fitness is the
% smaller the better.

[popSize,objDim]=size(objVals);
strength=zeros(popSize,1);
rawFit=zeros(popSize,1);
% Calculate the strength value---------------------------------------------
for i=1:popSize
    val=objVals(i,1:objDim);
    repObjv=val(ones(popSize,1),1:objDim);
    tmpRes=(objVals>=repObjv);% Temporary results
    tmpRes=sum(tmpRes,2);
    strength(i)=sum(tmpRes==objDim)-1;
end
% Calculate the raw fitness------------------------------------------------
for i=1:popSize
    val=objVals(i,1:objDim);
    repObjv=val(ones(popSize,1),1:objDim);
    tmpRes=(objVals<=repObjv);
    tmpRes=sum(tmpRes,2);
    rawFit(i)=sum(strength(tmpRes==objDim,1))-strength(i);
end
end