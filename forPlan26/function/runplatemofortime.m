%platemo('problem',@test1GOCEA_new,'algorithm',@GLT1,'D',3,'N',100,'maxFE',10000,'save',20);

% problem={@GLT1,@GLT2,@GLT3,@GLT4,@GLT5,@GLT6,@LZ1,@LZ2,@LZ3,@LZ4,@LZ5,@LZ6,@LZ7,@LZ8,@LZ9};
% problem={@planner_original};
problem={@GLT5};
probNum=size(problem,2);

% algorithm={@GOCEA,@NSGAIII,@PreDEMO,@PiMOEA,@MOEADSS,@IMMOEA,@RMMEDA,@SMEA};
%algorithm={@IMMOEA,@GOCEA,@NSGAII,@SMSEMOA,@MOEADCMA,@MOEADAWA,@RMMEDA,@SMEA,@OCEA};
%algorithm={@MOPSO,@rNSGAII,@WOF,@VaEA,@SPEA2,@SPEAR,@SMSEMOA,@SMEA,@MOEADD,@MOEADAWA,@MOEAD,@IMMOEAD,@HypE};
%algorithm={@MOEADSS};
algorithm={@SMEA};
algNum=size(algorithm,2);

vardim=[10];
%vardim=[30];
%M=[2,2,2,2,2,3,2,2,2];
% M=[2,2,2,2,3,3,2,2,2,2,2,3,2,2,2];
M=[3];
%N=[100,100,100,100,105,105,100,100,100,100,100,100,100,105,105,105,100,100,100,105,105,105,105,105,100,100,100,100,100];
%N=[100,100,100,100,100,105,100,100,100];
% N=[100,100,100,100,105,105,100,100,100,100,100,105,100,100,100];

N=[105];
runTimes=1;
workerNum=31;
maxgens=300;
% maxgens1=3000;
%maxFE=30000;

s=20;   %runplatemo函数中save值
A=zeros(1,algNum);
% B=zeros(9,algNum);
for m=1:probNum
    for n=1:algNum
        maxFE=maxgens*N(m);
        for q=1:31
        tic
        platemo('problem',problem{m},'algorithm',algorithm{n},'D',vardim,'N',N(m),'M',M(m),'maxFE',maxFE,'save',s);
        A(m,n)=toc+A(m,n);
        q
        end
    end
end
% for m=7:probNum
%     for n=1:algNum
%         maxFE=maxgens1*N(m);
%         for q=1:31
%         tic
%         platemo('problem',problem{m},'algorithm',algorithm{n},'D',vardim,'N',N(m),'M',M(m),'maxFE',maxFE,'save',s);
%         B(m-6,n)=toc+B(m-6,n);
%         q
%         end
%     end
% end
save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/DATAforSHENGAO/', ['_glt5smea','.mat']],'A')
% save(['/media/haichao/9c7b702c-d542-4ac0-8f33-8aac10340a9e/haichao/Desktop/wanghaoyue1/DATAforSHENGAO/', ['_LZ','.mat']],'B')
