
%problem={@WFG1,@WFG2,@WFG3,@WFG4,@WFG5,@WFG6,@WFG7,@WFG8,@WFG9};%可能是多个问题
%problem={@GLT1,@GLT2,@GLT3,@GLT4,@GLT5,@GLT6};
%problem={@LZ1,@LZ2,@LZ3,@LZ4,@LZ5,@LZ6,@LZ7,@LZ8,@LZ9};
problem={@planner_original};
probNum=size(problem,2);
N=[100];
algorithm={@GOCEA};%%%%%%%%%%
algNum=size(algorithm,2);
D=[30,60,90];
nD=size(D,2);
M=[3];

Kmax=[5,7,10,11,13];
nKmax=size(Kmax,2);

%BETA=[0.1,0.3,0.5,0.7,0.9];
BETA=[0.8,0.9];
nBETA=size(BETA,2);

%F=[0.1,0.3,0.5,0.7,0.9];
F=[0.5,0.7];
nF=size(F,2);

%CR=[0.2,0.4,0.6,0.8,1];
CR=[0.6,0.8,1];
nCR=size(CR,2);

runTimes=1;
workerNum=31;
maxgens=300;
%maxFE=30000;
s=20;
for m=1:probNum
for d=1:nD
    for i=1:nKmax
        for k=1:nF
            for l=1:nCR
                for p=1:nBETA
                    for n=1:algNum
                    maxFE=maxgens*N;
                    GOCEASpmdRun(problem{m},algorithm{n},Kmax(i),BETA(p),F(k),CR(l),N,M,D(d),maxFE,s,runTimes,workerNum);
                    end
                end
            end
        end
    end
end
end
