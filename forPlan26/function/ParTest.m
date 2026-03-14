% function ParTest_GLTa
problem={@LZ6,@LZ7,@LZ8,@LZ9};%可能是多个问题
probNum=size(problem,2);
N=[105,100,100,100];
algorithm={@GOCEA,@NSGAII,@SMSEMOA,@IMMOEA,@MOEAD,@MOEADCMA,@MOEADAWA,@RMMEDA,@SMEA};
algNum=size(algorithm,2);
D=[30];
nD=size(D,2);
M=[3,2,2,2];
runTimes=1;
workerNum=31;
maxFE=60000;
s=20;
for d=1:nD
    for m=1:probNum
        for n=1:algNum
        SpmdRun(problem{m},algorithm{n},N(m),M(m),D(d),maxFE,s,runTimes,workerNum);
        end
    end
end
% end