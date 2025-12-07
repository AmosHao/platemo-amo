function newPop=testPolynomialMutation(pop,bounds,pm)
% This function implements a polynomial mutation operator. The codes have
% been checked by Hu Zhang on April 15th 2015.
% Author: Zhang Hu, Harbin Institute of Technology
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)
%%
    [popSize,varDim]=size(pop);
    newPop=zeros(popSize,varDim);
    for i=1:popSize
        sol=pop(i,1:varDim);
        newPop(i,1:varDim)=Mutation(sol,bounds,pm);
    end
end

function newSol=Mutation(sol,bounds,pm)
    eta_m=20; % distribution index
    varDim=size(sol,2);
    newSol=sol;
    for i=1:varDim
        if rand<=pm
            y=sol(1,i);
            yL=bounds(1,i);
            yu=bounds(2,i);
            delta1=(y-yL)/(yu-yL);
            delta2=(yu-y)/(yu-yL);
            rnd=rand;
            mut_pow=1.0/(eta_m+1.0);
            if rnd<=0.5
                xy=1.0-delta1;
                val= 2.0*rnd+(1.0-2.0*rnd)*(xy^(eta_m+1.0));
                deltaq=(val^mut_pow)-1.0;
            else
                xy=1.0-delta2;
                val= 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(xy^(eta_m+1.0));
                deltaq=1.0-(val^mut_pow);
            end
            y=y+deltaq*(yu-yL);
            if (y<yL)
                % y=yL-rand*(y-yL);
                y=yL;
            end
            if (y>yu)
                % y=yu+rand*(yu-y);
                y=yu;
            end
            newSol(1,i)=y;
        end
    end
end