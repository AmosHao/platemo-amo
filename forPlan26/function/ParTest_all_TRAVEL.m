% function ImParTest
%problem={@GLT1,@GLT2,@GLT3,@GLT4,@GLT5,@GLT6,@UF1,@UF2,@UF3,@UF4,@UF5,@UF6,@UF7,@UF8,@UF9,@UF10,@IMOP1,@IMOP2,@IMOP3,@IMOP4,@IMOP5,@IMOP6,@IMOP7,@IMOP8,@ZDT1,@ZDT2,@ZDT3,@ZDT4,@ZDT6};
%problem={@WFG1,@WFG2,@WFG3,@WFG4,@WFG5,@WFG6,@WFG7,@WFG8,@WFG9};
%problem={@LZ1,@LZ2,@LZ3,@LZ4,@LZ5,@LZ6,@LZ7,@LZ8,@LZ9};
%problem={@GLT3,@GLT1,@GLT2,@GLT4,@GLT5,@GLT6};
problem={@planner_points};
probNum=size(problem,2);

% algorithm={@GOCEA,@NSGAIII,@PreDEMO,@PiMOEA,@PeEA,@MOEADSS,@RMMEDA,@IMMOEA,@SMEA};
%algorithm={@IMMOEA,@GOCEA,@NSGAII,@SMSEMOA,@MOEADCMA,@MOEADAWA,@RMMEDA,@SMEA,@OCEA};
%algorithm={@MOPSO,@rNSGAII,@WOF,@VaEA,@SPEA2,@SPEAR,@SMSEMOA,@SMEA,@MOEADD,@MOEADAWA,@MOEAD,@IMMOEAD,@HypE};
algorithm={@GOCEAfortravel_singlepath_3};
%algorithm={@RMMEDA,@IMMOEA};
algNum=size(algorithm,2);

% vardim=[30];
%vardim=[30];
%M=[2,2,2,2,2,3,2,2,2];
%M=[2,2,2,2,3,3];
M=[3];
%N=[100,100,100,100,105,105,100,100,100,100,100,100,100,105,105,105,100,100,100,105,105,105,105,105,100,100,100,100,100];
%N=[100,100,100,100,100,105,100,100,100];
%N=[100,100,100,100,105,105];

N=[100];
runTimes=1;
workerNum=11;
% maxgens=300;
maxFE=20000;

s=20;   %runplatemo函数中save值

dotss=[0.2,0.2,0.8;
    2.7,0.5,0.8;
    3.5,3,0.8;
    2.9,4.8,0.8;
    2.5,6.8,0.8;
    0.5,5,0.8];
combines=nchoosek(1:size(dotss,1),2);
combinenum=size(combines,1);
for m=1:probNum
    for n=1:algNum
        for k=1:combinenum
            dots=[dotss(combines(k,1),:),dotss(combines(k,2),:)];
            dot1=combines(k,1);
            dot2=combines(k,2);
        % maxFE=maxgens*N(m);
        % disp(problem{m},algorithm{n},N(m),maxFE,s,runTimes,workerNum,M(m),parameter,dot1,dot2);
         % SpmdRun_all_TRAVEL(problem{m},algorithm{n},N(m),maxFE,s,runTimes,workerNum,M(m),parameter,dot1,dot2);
         SpmdRun_all_TRAVEL(problem{m},algorithm{n},N(m),maxFE,s,runTimes,workerNum,M(m),dots,dot1,dot2);
        end
    end
end
              