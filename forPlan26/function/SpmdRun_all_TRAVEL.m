function SpmdRun_all_TRAVEL(problem,algorithm,N,maxFE,s,runTimes,workerNum,M,dots,dot1,dot2)

% This script demonstrates using one EMO algorithm to compute
% multiple problems.
% runTimes=3;                     % The run times of each lab
% workerNum=10;                    % The number of labs
spmd
 
    p = struct();% 创建一个空的结构体
    q = struct();
    w = struct();

    q.Obj=cell(runTimes,1);
    q.Dec=cell(runTimes,1);
    q.Dec_end=cell(runTimes,1);
    q.validtrait_end=cell(runTimes,1);
    q.validtrait=cell(runTimes,1);

    p.runplatemo.igdv = cell(runTimes, 1);
    p.runplatemo.ihv = cell(runTimes, 1);
    p.runplatemo.Obj_end = cell(runTimes, 1);
    p.runplatemo.igdv_end = cell(runTimes, 1);
    p.runplatemo.ihv_end = cell(runTimes, 1);

    w.clustName=cell(runTimes,1);
    w.clustTag=cell(runTimes,1);
    w.clustName_end=cell(runTimes,1);
    w.clustTag_end=cell(runTimes,1);

    for j=1:runTimes
        
       [validtrait_end,validtrait,clustName,clustTag,clustName_end,clustTag_end,Dec,Obj,Dec_end,Obj_end,igd,hv,igd_end,hv_end]=runplatemo_all_points(problem,algorithm,N,maxFE,s,M,dots);
      
        p.runplatemo.igdv(j,1)={igd};
        p.runplatemo.ihv(j,1)={hv};
        p.runplatemo.Obj_end(j,1)={Obj_end};
        p.runplatemo.igdv_end(j,1)={igd_end};
        p.runplatemo.ihv_end(j,1)={hv_end};
        % p.runplatemo.Objgens(j,1)={objgens};
        % p.runplatemo.Popgens(j,1)={popgens};
        
        q.Obj{j,1}={Obj};
        q.Dec{j,1}={Dec};
        q.Dec_end{j,1}={Dec_end};
        q.validtrait_end{j,1}={validtrait_end};
        q.validtrait{j,1}={validtrait};

        w.clustName{j,1}={clustName};
        w.clustTag{j,1}={clustTag};
        w.clustName_end{j,1}={clustName_end};
        w.clustTag_end{j,1}={clustTag_end};

    end
end
% Assemble the results from the labs
fun=p{1};
qun=q{1};
wun=w{1};
for i=2:workerNum
    tmpFun=p{i};
    tmpQun=q{i};
    tmpWun=w{i};

    fun.runplatemo.igdv=[fun.runplatemo.igdv;tmpFun.runplatemo.igdv];
    fun.runplatemo.ihv=[fun.runplatemo.ihv;tmpFun.runplatemo.ihv];
    fun.runplatemo.Obj_end=[fun.runplatemo.Obj_end;tmpFun.runplatemo.Obj_end];
    fun.runplatemo.igdv_end=[fun.runplatemo.igdv_end;tmpFun.runplatemo.igdv_end];
    fun.runplatemo.ihv_end=[fun.runplatemo.ihv_end;tmpFun.runplatemo.ihv_end];

     qun.Obj=[qun.Obj;tmpQun.Obj];
     qun.Dec=[qun.Dec;tmpQun.Dec];
     qun.Dec_end=[qun.Dec_end;tmpQun.Dec_end];
     qun.validtrait_end=[qun.validtrait_end;tmpQun.validtrait_end];
     qun.validtrait=[qun.validtrait;tmpQun.validtrait];

     wun.clustName=[wun.clustName;tmpWun.clustName];
     wun.clustTag=[wun.clustTag;tmpWun.clustTag];
     wun.clustName_end=[wun.clustName_end;tmpWun.clustName_end];
     wun.clustTag_end=[wun.clustTag_end;tmpWun.clustTag_end];
end
prob=fun;
qrob=qun;
wrob=wun;


%maxGens=num2str(maxGens);

info=functions(problem);
probleminfo=info.function;
info=functions(algorithm);
algorithminfo=info.function;

N=num2str(N);
% vardim=num2str(vardim);
maxFE=num2str(maxFE);
runT=num2str(runTimes*workerNum);
dot1=num2str(dot1);
dot2=num2str(dot2);

save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/datafortravel/hangji_newdots/1/', ['_problem',probleminfo,'_algorithm',algorithminfo,'_N',N,'dot1',dot1,'dot2',dot2,'_maxFE',maxFE,'_runT',runT,'.mat']],'prob');
save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/datafortravel/hangji_newdots/2/', ['_problem',probleminfo,'_algorithm',algorithminfo,'_N',N,'dot1',dot1,'dot2',dot2,'_maxFE',maxFE,'_runT',runT,'.mat']],'qrob');
save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/datafortravel/hangji_newdots/3/', ['_problem',probleminfo,'_algorithm',algorithminfo,'_N',N,'dot1',dot1,'dot2',dot2,'_maxFE',maxFE,'_runT',runT,'.mat']],'wrob');
%save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/PlatEMO-master-using/rundatanew/testOCEAonGLT/k7_1/', ['_problem',probleminfo,'_algorithm',algorithminfo,'_N',N,'_D',vardim,'_maxFE',maxFE,'_runT',runT,'.mat']],'prob');
%save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/PlatEMO-master-using/rundatanew/testOCEAonGLT/k7_2/', ['_problem',probleminfo,'_algorithm',algorithminfo,'_N',N,'_D',vardim,'_maxFE',maxFE,'_runT',runT,'.mat']],'qrob');
%save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/PlatEMO-master-using/rundatanew/testOCEAonGLT/k7_3/', ['_problem',probleminfo,'_algorithm',algorithminfo,'_N',N,'_D',vardim,'_maxFE',maxFE,'_runT',runT,'.mat']],'wrob');
%save(['C:\Users\DELL\Desktop\PlatEMO-master\PlatEMO-master\PlatEMO\mydate\partest', ['result','_problem',probleminfo,'_algorithm',algorithminfo,'_Kmax',Kmax,'_BETA',BETA,'_F',F,'_CR',CR,'_run',runT,'_maFE',maxFE,'.mat']],'result');
clear fun
end


