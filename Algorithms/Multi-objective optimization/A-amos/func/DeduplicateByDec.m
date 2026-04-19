function Population = DeduplicateByDec(Population)
%DeduplicateByDec - Remove duplicate solutions by decision variables.
%
%   Population = DeduplicateByDec(Population) keeps only one solution for
%   each unique decision vector (Population.decs). For duplicates, it keeps
%   the solution with the smallest objective 1; ties are broken by objective
%   2, then objective 3, etc. (lexicographic on objectives).
%
%   Note: This function does NOT evaluate solutions. It assumes .decs and
%   .objs are already available (SOLUTION/INDIVIDUAL array).

    if isempty(Population)
        return;
    end
    Decs = Population.decs;
    Objs = Population.objs;
    if isempty(Decs) || isempty(Objs)
        return;
    end

    % Group by identical decision vectors (stable order for reproducibility)
    [~,~,gid] = unique(Decs,'rows','stable');
    K = max(gid);
    keep = false(size(gid));

    M = size(Objs,2);
    sortCols = 1:M;
    for g = 1:K
        idx = find(gid==g);
        if numel(idx) == 1
            keep(idx) = true;
        else
            [~,ord] = sortrows(Objs(idx,:), sortCols);
            keep(idx(ord(1))) = true;
        end
    end
    Population = Population(keep);
end

