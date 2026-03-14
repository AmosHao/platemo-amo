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
workerNum=31;
% maxgens=300;
maxFE=30000;

s=20;   %runplatemo函数中save值


for m=1:probNum
    for n=1:algNum
        maxFE=maxgens*N(m);
         SpmdRun_all(problem{m},algorithm{n},N(m),maxFE,s,runTimes,workerNum,vardim,M(m));
    end
end
              