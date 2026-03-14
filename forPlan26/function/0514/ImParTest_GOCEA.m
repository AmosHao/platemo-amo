% function ImParTest
%problem={@GLT1,@GLT2,@GLT3,@GLT4,@GLT5,@GLT6};
problem={@LZ7,@LZ6,@LZ8,@LZ9,@LZ1,@LZ2,@LZ3,@LZ4,@LZ5};
probNum=size(problem,2);
% N=[100,100,100,100,105,105];
N=[100,300,100,100,100,100,100,100,100];
algorithm={@GOCEA};%%%%%%%%%%
algNum=size(algorithm,2);
vardim=[30];
nD=size(vardim,2);
%M=[2,2,2,2,3,3];
M=[2,3,2,2,2,2,2,2,2];
Kmax=[3,5,7,9,11];
% Kmax=[7];
nKmax=size(Kmax,2);

%BETA=[0.1,0.3,0.5,0.7,0.9];
BETA=[0.9];
nBETA=size(BETA,2);

%F=[0.1,0.3,0.5,0.7,0.9];
F=[0.7];
nF=size(F,2);

%CR=[0.2,0.4,0.6,0.8,1];
CR=[1];
nCR=size(CR,2);

runTimes=1;
workerNum=31;
maxFE=[300000,1000000,300000,300000,300000,300000,300000,300000,300000];
s=20;
% for i=1:probNum
%     if objDim(i)==2
%         N(i)=100;
%     else
%         N(i)=105;
%     end
% end
for m=1:probNum
for w=1:nKmax
for i=1:nBETA
    for k=1:nF
        for l=1:nCR
                    for n=1:algNum
                    ImSpmdRun_GOCEA(problem{m},algorithm{n},Kmax(w),BETA(i),F(k),CR(l),N(m),maxFE(m),s,runTimes,workerNum,vardim,M(m));
                    end
                end
               
        end
    end
end
end