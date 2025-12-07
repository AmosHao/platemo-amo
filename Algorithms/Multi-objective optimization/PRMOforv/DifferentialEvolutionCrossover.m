function offs = DifferentialEvolutionCrossover(parents,bounds,F,CR)
% Author:  Hu Zhang, Harbin Institute of Technology
% Last modified: October 26, 2013
% Copyright (C) 2013-2016 by Hu Zhang (e-mail: jxzhanghu@126.com)

if nargin==1
    error('input is not enough!');
elseif nargin==2
    F=0.5;
    CR=1;
end
varDim=size(parents,2);
% varDim=3*n;
offs=parents(3,1:varDim);
%     randInt=round(1+(variableCount-1)*rand(1));
for i=1:varDim
    rnd = rand;
    %         if (rnd<CR)||(i==randInt)
    if (rnd<CR)
        value = parents(3,i)+F*(parents(1,i)-parents(2,i));
        if value<bounds(1,i)
           value=bounds(1,i);
            %                 value=bounds(1,i)+rand*(bounds(2,i)-bounds(1,i));
        elseif value>bounds(2,i)
            value=bounds(2,i);
            %                 value=bounds(2,i)-rand*(bounds(2,i)-bounds(1,i));
        end
        offs(1,i)=value;
    end
end
end