% function ImParTest
%problem={@GLT1,@GLT2,@GLT3,@GLT4,@GLT5,@GLT6,@UF1,@UF2,@UF3,@UF4,@UF5,@UF6,@UF7,@UF8,@UF9,@UF10,@IMOP1,@IMOP2,@IMOP3,@IMOP4,@IMOP5,@IMOP6,@IMOP7,@IMOP8,@ZDT1,@ZDT2,@ZDT3,@ZDT4,@ZDT6};
%problem={@WFG1,@WFG2,@WFG3,@WFG4,@WFG5,@WFG6,@WFG7,@WFG8,@WFG9};
%problem={@LZ1,@LZ2,@LZ3,@LZ4,@LZ5,@LZ6,@LZ7,@LZ8,@LZ9};
%problem={@GLT3,@GLT1,@GLT2,@GLT4,@GLT5,@GLT6};
problem={@planner_travel_maxhjd_newobj};
% problem={@planner_simple_checkrestrain_2_con};
probNum=size(problem,2);

% algorithm={@GOCEA,@NSGAIII,@PreDEMO,@PiMOEA,@PeEA,@MOEADSS,@RMMEDA,@IMMOEA,@SMEA};
%algorithm={@IMMOEA,@GOCEA,@NSGAII,@SMSEMOA,@MOEADCMA,@MOEADAWA,@RMMEDA,@SMEA,@OCEA};
%algorithm={@MOPSO,@rNSGAII,@WOF,@VaEA,@SPEA2,@SPEAR,@SMSEMOA,@SMEA,@MOEADD,@MOEADAWA,@MOEAD,@IMMOEAD,@HypE};
% algorithm={@MOEADSS,@NSGAIII,@PiMOEA,@PeEA,@GOCEA};
% algorithm={@PRMOEADSS,@PRNSGAIII,@PRPiMOEA,@PRPeEA,@PRGOCEA};
% algorithm={@CMOEA_MS,@CMOSMA,@CTSEA};
% algorithm={@PRCMOEA_MS,@PRCMOSMA,@PRCTSEA};
% algorithm={@PRMO,@PRMO_generate,@PRMO_beta,@PRMO_repair};
% algorithm={@simplePRMO};
% algorithm={@simplePRMOEADSS,@simplePRNSGAIII,@simplePRPiMOEA,@simplePRPeEA,@simplePRGOCEA};
% algorithm={@simplePRCMOEA_MS,@simplePRCMOSMA};
%algorithm={@RMMEDA,@IMMOEA};
% algorithm={@PRPiMOEA,@PRNSGAIII,@PRCMOEA_MS,@PRCMOSMA,@PRPeEA,@PRMOEADSS,@PRGOCEA};
% algorithm={@PRMO_OCmove,@PRMO};
algorithm={@RMMEDA,@RVEA,@LMEA};
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

% dotss=[0.2,0.2,0.8;
%     2.7,0.5,0.8;
%     3.5,3,0.8;
%     2.9,4.8,0.8;
%     2.5,6.8,0.8;
%     0.5,5,0.8];
% combines=nchoosek(1:size(dotss,1),2);
% combinenum=size(combines,1);
time=zeros(probNum,algNum);
for m=1:probNum
    tic;
    for n=1:algNum
        % for
        % k=1:combinenum
            % dots=[dotss(combines(k,1),:),dotss(combines(k,2),:)];
            dots=[500,500,800,7000,3000,800];
            % dots=[500,500,800,600,600,800];
            % dots=[50,50,800,600,600,800];
        %     dot1=combines(k,1);
        %     dot2=combines(k,2);
        % maxFE=maxgens*N(m);
        % disp(problem{m},algorithm{n},N(m),maxFE,s,runTimes,workerNum,M(m),parameter,dot1,dot2);
         % SpmdRun_all_TRAVEL(problem{m},algorithm{n},N(m),maxFE,s,runTimes,workerNum,M(m),parameter,dot1,dot2);
         SpmdRun_all_TRAVEL_PRMO_class(problem{m},algorithm{n},N(m),maxFE,s,runTimes,workerNum,M(m),dots);
        % end
    atime=toc;
    time(m,n)=atime;
    end

end

save(['/home/haichao/Documents/why/wanghaoyue1/20260314OUPAcompare', ['_OCmovePRMOtime','.mat']],'time');