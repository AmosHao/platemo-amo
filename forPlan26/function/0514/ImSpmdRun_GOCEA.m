function ImSpmdRun_GOCEA(problem,algorithm,K,BETA,F,CR,N,maxFE,s,runTimes,workerNum,vardim,M)

% This script demonstrates using one EMO algorithm to compute
% multiple problems.
% runTimes=3;                     % The run times of each lab
% workerNum=10;                    % The number of labs
spmd
 
    p = struct();% 创建一个空的结构体

    p.runplatemo.igdv = cell(runTimes, 1);
    p.runplatemo.ihv = cell(runTimes, 1);
    p.runplatemo.Obj = cell(runTimes, 1);
    p.runplatemo.igdv_end = cell(runTimes, 1);
    p.runplatemo.ihv_end = cell(runTimes, 1);
    for j=1:runTimes
        
       [Obj,igd,hv,igd_end,hv_end]=runplatemo_GOCEA_OLD(problem,algorithm,K,BETA,F,CR,N,maxFE,s,vardim,M);
      
        p.runplatemo.igdv(j,1)={igd};
        p.runplatemo.ihv(j,1)={hv};
        p.runplatemo.Obj(j,1)={Obj};
        p.runplatemo.igdv_end(j,1)={igd_end};
        p.runplatemo.ihv_end(j,1)={hv_end};

    end
end
% Assemble the results from the labs
fun=p{1};
for i=2:workerNum
    tmpFun=p{i};
    fun.runplatemo.igdv=[fun.runplatemo.igdv;tmpFun.runplatemo.igdv];
    fun.runplatemo.ihv=[fun.runplatemo.ihv;tmpFun.runplatemo.ihv];
    fun.runplatemo.Obj=[fun.runplatemo.Obj;tmpFun.runplatemo.Obj];
    fun.runplatemo.igdv_end=[fun.runplatemo.igdv_end;tmpFun.runplatemo.igdv_end];
    fun.runplatemo.ihv_end=[fun.runplatemo.ihv_end;tmpFun.runplatemo.ihv_end];
end
prob=fun;

runT=num2str(runTimes*workerNum);
%maxGens=num2str(maxGens);
F=num2str(F);
CR=num2str(CR);
info=functions(problem);
probleminfo=info.function;
info=functions(algorithm);
algorithminfo=info.function;
BETA=num2str(BETA);
K=num2str(K);
N=num2str(N);
M=num2str(M);
vardim=num2str(vardim);
maxFE=num2str(maxFE);


save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/PlatEMO-master-using/yuyu/MYDATA/testLZpra0514/', ['_problem',probleminfo,'_algorithm',algorithminfo,'_Kmax',K,'_BETA',BETA,'_F',F,'_CR',CR,'_N',N,'M',M,'_D',vardim,'_maxFE',maxFE,'_runT',runT,'.mat']],'prob');
%save(['C:\Users\DELL\Desktop\PlatEMO-master\PlatEMO-master\PlatEMO\mydate\partest', ['result','_problem',probleminfo,'_algorithm',algorithminfo,'_Kmax',Kmax,'_BETA',BETA,'_F',F,'_CR',CR,'_run',runT,'_maFE',maxFE,'.mat']],'result');
clear fun
end


