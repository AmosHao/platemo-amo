function offs=SBXCrossover(parents,bounds,pc)
% This function allows to apply a simulated binary crossover (SBX) operator
% using two parent solutions. These codes have been checked by Hu Zhang on
% April 15th 2015.
% Author: Zhang Hu, Harbin Institute of Technology
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

    %%
    eta_c=20;               % Distribution index
    epsilon=1.0e-14;        % The minimum difference between two real values
    if nargin<3
        pc=0.9;             % The default Pc is 0.9 
    end
    varDim=size(parents,2); % The number of decision variables
    offs=parents;
    r1=rand;
    if r1<=pc 
        for i=1:varDim
            x1=parents(1,i);x2=parents(2,i);
            r2=rand;
            if r2<=0.5
                if abs(x1-x2)>epsilon
                    if x1<x2
                        y1=x1;y2=x2;
                    else
                        y1=x2;y2=x1;
                    end
                    yL=bounds(1,i);yu=bounds(2,i);
                    r3=rand;
                    beta=1.0+(2.0*(y1-yL)/(y2-y1));alpha=2.0-beta^(-(eta_c+1.0));
                    if r3<=(1.0/alpha)
                        betaq=(r3*alpha)^(1.0/(eta_c+1.0));
                    else
                        betaq=(1.0/(2.0-r3*alpha))^(1.0/(eta_c+1.0));
                    end
                    c1=0.5*((y1+y2)-betaq*(y2-y1));
                    beta=1.0+(2.0*(yu-y2)/(y2-y1));alpha=2.0-(beta^(-(eta_c+1.0)));
                    if r3<=(1.0/alpha)
                        betaq=(r3*alpha)^(1.0/(eta_c+1.0));
                    else
                        betaq=(1.0/(2.0-r3*alpha))^(1.0/(eta_c+1.0));
                    end
                    c2=0.5*((y1+y2)+betaq*(y2-y1));
                    c1=max(c1,yL);c2=max(c2,yL);c1=min(c1,yu);c2=min(c2,yu);
                    if rand<=0.5
                        offs(1,i)=c2;offs(2,i)=c1;
                    else
                        offs(1,i)=c1;offs(2,i)=c2;
                    end 
                else
                    offs(1,i)=x1;offs(2,i)=x2;
                end
            else
                offs(1,i) = x2;offs(2,i) = x1;
            end
        end
    end
end