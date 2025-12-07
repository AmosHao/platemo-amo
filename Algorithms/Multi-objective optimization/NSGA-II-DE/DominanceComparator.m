function d = DominanceComparator(fA, fB)
% Compares two solutions A and B given their objective function
% values fA and fB. Returns whether A dominates B.
%
% Input:
% - fA					- The objective function values of solution A
% - fB					- The objective function values of solution B
%
% Output:
% - d					- d is 1 if fA dominates fB; d is -1 if fB
%                         dominates fA, otherwise d equals to 0.
%
% Author: Zhang Hu, Harbin Institute of Technology
% Last modified: October 27, 2013
% Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

    %%
    % Elegant, but not very efficient
    d=0;
    if (all(fA <= fB) && any(fA < fB))
        d = 1;
    elseif (all(fB <= fA) && any(fB < fA))
        d = -1;
    end
end