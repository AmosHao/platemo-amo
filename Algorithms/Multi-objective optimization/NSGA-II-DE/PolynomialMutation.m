function offs = PolynomialMutation(parents,bounds,Pm)
% This function implements a polynomial mutation operator.
%
% Input:
% - parent		    - The solution to join the mutation
% - bounds          - The bounds of the decision variables
% - Pm              - The mutation probability
%
% Output:
% - offs			- The solution ofter mutation
%
% Author: Zhang Hu, Harbin Institute of Technology
% Last modified: October 27, 2013
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

    [pSize,pCount]=size(parents);
    offs=zeros(pSize,pCount);
    for i=1:pSize
        p=parents(i,1:pCount);
        offs(i,1:pCount)=Mutation(p,bounds,Pm);
    end
end

function offs = Mutation(parent,bounds,Pm)
    eta=20; % distribution index
    variableCount=size(parent,2);
    offs=parent;

    for i=1:variableCount
        if rand<= Pm
            y=parent(1,i);
            yL=bounds(1,i);
            yu=bounds(2,i);
            delta1=(y-yL)/(yu-yL);
            delta2=(yu-y)/(yu-yL);
            rnd=rand;
            mut_pow=1.0/(eta+1.0);
            if rnd<=0.5
                xy=1.0-delta1;
                val= 2.0*rnd+(1.0-2.0*rnd)*(xy^(eta+1.0));
                deltaq=(val^mut_pow)-1.0;
            else
                xy=1.0-delta2;
                val= 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(xy^(eta+1.0));
                deltaq=1.0-(val^mut_pow);
            end
            y=y+deltaq*(yu-yL);
            if (y<yL)
                y=yL-rand*(y-yL);
            end
            if (y>yu)
                y=yu+rand*(yu-y);
            end
            offs(1,i)=y;
        end
    end
end