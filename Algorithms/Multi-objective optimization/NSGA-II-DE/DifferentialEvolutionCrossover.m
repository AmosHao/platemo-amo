function childIndividual = DifferentialEvolutionCrossover(parents,bounds,F,CR)
% This function implements the DE crossover operator
%
% Input:
% - parents      - the solutions used to join DE operation, it has
% three rows, and each row represents a solution
% - bounds       - the bounds of the variables, the 1st row is the
% lower bound, and the 2nd row is the upper bound
% - F,CR         - the control parameters in DE
%
% Output:
% -childIndividual - the generated individual after DE crossover operation
%
% Author: Zhang Hu, Harbin Institute of Technology
% Last modified: October 26, 2013
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

%%
    if nargin==1
        error('input is not enough!');
    elseif nargin==2
        F=0.5;
        CR=1;
    end
    variableCount=size(parents,2);
    % The 3rd individual included in parents is the current individual
    childIndividual=parents(3,1:variableCount);
    randInt=randsample(variableCount,1);
    rnd = rand(1,variableCount);
    for i=1:variableCount
        if (rnd(i)<CR)||(i==randInt)
            value = parents(3,i)+F*(parents(1,i)-parents(2,i));
            % if value<bounds(1,i)
            %     value=bounds(1,i);
            %     %                 value=bounds(1,i)-rand*(value-bounds(1,i));
            % elseif value>bounds(2,i)
            %     value=bounds(2,i);
            %     %                 value=bounds(2,i)+rand*(bounds(2,i)-value);
            % end
            childIndividual(1,i)=value;
        end
    end
end