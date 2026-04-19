function Dec = RepairOrderEncoding(Problem,Dec)
%RepairOrderEncoding - Minimal repair for order-assignment sequence encoding.
%
% This is a baseline feasibility adapter for problems whose decision vector
% represents:
%   - points 1..2*n (each should appear once)
%   - (m-1) zeros as separators
% and each merchant i must appear before its paired customer (i+n).
%
% The repair is intentionally generic/minimal: it ensures integrality,
% bounds, completeness, correct number of zeros, non-empty segments (best
% effort), and preserves as much ordering signal as possible.

    if isempty(Dec)
        return;
    end
    Dec = round(Dec);

    n = Problem.n;
    m = Problem.m;
    D = Problem.D;
    maxPoint = 2*n;

    % Clamp to feasible numeric range
    Dec = min(max(Dec,0),maxPoint);

    N = size(Dec,1);
    for r = 1:N
        row = Dec(r,:);

        % Derive a stable merchant order from the current row (best effort)
        merchants = row(row>=1 & row<=n);
        merchants = merchants(:)';
        merchants = unique(merchants,'stable');
        missingM = setdiff(1:n, merchants, 'stable');
        if ~isempty(missingM)
            merchants = [merchants missingM(randperm(numel(missingM)))];
        end
        if numel(merchants) > n
            merchants = merchants(1:n);
        end

        % Build paired sequence [i, i+n] to satisfy precedence
        seq = zeros(1,2*n);
        for k = 1:n
            mi = merchants(k);
            seq(2*k-1) = mi;
            seq(2*k)   = mi + n;
        end

        % Determine separator insertion positions (must be after even indices)
        targetZeros = max(0,m-1);
        if targetZeros == 0
            repaired = seq;
        else
            evenPositions = 2:2:(2*n-2); % valid cut points
            zpos = find(row==0);
            % Map existing zero positions to nearest valid even cut points
            if ~isempty(zpos) && ~isempty(evenPositions)
                mapped = zeros(1,numel(zpos));
                for i = 1:numel(zpos)
                    [~,ii] = min(abs(evenPositions - zpos(i)));
                    mapped(i) = evenPositions(ii);
                end
                mapped = unique(mapped,'stable');
            else
                mapped = [];
            end

            % Ensure exactly targetZeros cut points
            mapped = mapped(mapped>=2 & mapped<=2*n-2);
            mapped = unique(mapped,'stable');
            if numel(mapped) > targetZeros
                mapped = mapped(1:targetZeros);
            elseif numel(mapped) < targetZeros
                pool = setdiff(evenPositions, mapped, 'stable');
                if ~isempty(pool)
                    add = pool(randperm(numel(pool), min(numel(pool), targetZeros-numel(mapped))));
                    mapped = sort([mapped add]);
                end
            end

            % Insert zeros after selected cut points
            repaired = [];
            prev = 0;
            mapped = sort(mapped);
            for i = 1:numel(mapped)
                repaired = [repaired seq(prev+1:mapped(i)) 0]; %#ok<AGROW>
                prev = mapped(i);
            end
            repaired = [repaired seq(prev+1:end)];
        end

        % Final length adjust (pad/truncate) to D
        if numel(repaired) < D
            repaired = [repaired zeros(1,D-numel(repaired))];
        elseif numel(repaired) > D
            repaired = repaired(1:D);
        end

        Dec(r,:) = repaired;
    end
end

