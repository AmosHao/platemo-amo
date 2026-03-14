
function SpmdRun(problem,algorithm,N,M,D,maxFE,s,runTimes,workerNum)
% This script demonstrates using one EMO algorithm to compute
% multiple problems.
% runTimes=3;                     % The run times of each lab
% workerNum=7;                    % The number of labs
spmd
    %p=MOPs(probName,varDim,objDim); % Determine the problem to be solved
    % 创建一个空的结构体
    p = struct();

% 为结构体添加字段
    p.runplatemo.igdv = cell(runTimes, 1);
    p.runplatemo.ihv = cell(runTimes, 1);
    p.runplatemo.Obj = cell(runTimes, 1);
    p.runplatemo.igdv_end = cell(runTimes, 1);
    p.runplatemo.ihv_end = cell(runTimes, 1);
    for j=1:runTimes
        % Generate objective values using APMOHV
        % addpath('E:\PlatEMO-master\PlatEMO\runplatemo.m');
       [Obj,igd,hv,igd_end,hv_end]=runplatemo(problem,algorithm,N,M,D,maxFE,s);
        % rmpath('E:\PlatEMO-master\PlatEMO\runplatemo.m');
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
maxFE=num2str(maxFE);
N=num2str(N);
M=num2str(M);
D=num2str(D);
info=functions(problem);
probleminfo=info.function;
info=functions(algorithm);
algorithminfo=info.function;


save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/PlatEMO-master-using/PlatEMO/rundata/GOCEA_experimence/GOCEA_testalgo_onLZ/', ['_problem',probleminfo,'_algorithm',algorithminfo,'_N',N,'_D',D,'M',M,'_maxFE',maxFE,'_runT',runT,'.mat']],'prob');
%save(['E:\PlatEMO-master\PlatEMO\testdataparallel\', ['result','_problem',probleminfo,'_algorithm',algorithminfo,'_Kmax',Kmax,'_BETA',BETA,'_F',F,'_CR',CR,'_run',runT,'_maFE',maxFE,'.mat']],'result');
clear fun
end


