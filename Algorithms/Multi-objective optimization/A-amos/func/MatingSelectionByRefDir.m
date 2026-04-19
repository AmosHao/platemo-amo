function [p1,p2,state] = MatingSelectionByRefDir(Population,N,Z,Zmin,varargin)
%MatingSelectionByRefDir - Mating selection by NSGA-III reference directions.
%
%   [p1,p2,state] = MatingSelectionByRefDir(Population,N,Z,Zmin,...) returns
%   two index vectors p1 and p2 (both length N) such that p2(i) is selected
%   from a *different* reference-direction niche than p1(i) whenever possible.
%
%   Inputs:
%     Population : INDIVIDUAL array
%     N          : number of mating pairs to output
%     Z          : reference points (NZ x M), from UniformPoint
%     Zmin       : ideal point (1 x M) used by NSGA-III normalization
%
%   Name-value options:
%     'p1'   : predefined p1 indices (default: 1:min(N,|Pop|), then random)
%     'mode' : 'random' | 'farthest' (default: 'farthest')
%              - 'farthest': choose p2 maximizing normalized objective distance
%     'state': struct for caching (currently stores pi/d/PopObjN)
%
%   Output state fields:
%     state.pi     : niche index (1..NZ) for each individual in Population
%     state.PopObj : normalized objectives used for association (|Pop| x M)
%
%   Notes:
%     - Uses the same normalization + association logic as NSGA-III
%       EnvironmentalSelection (LastSelection).
%

    popSize = numel(Population);
    if popSize == 0
        p1 = [];
        p2 = [];
        state = struct('pi',[],'PopObj',[]);
        return;
    end
    if nargin < 2 || isempty(N)
        N = popSize;
    end

    % ---- parse options ----
    p1 = [];
    mode = 'farthest';
    state = struct();
    for i = 1 : 2 : length(varargin)
        switch lower(varargin{i})
            case 'p1'
                p1 = varargin{i+1};
            case 'mode'
                mode = lower(string(varargin{i+1}));
            case 'state'
                state = varargin{i+1};
        end
    end

    % ---- build p1 ----
    if isempty(p1)
        if N <= popSize
            p1 = 1:N;
        else
            p1 = [1:popSize, randi(popSize,1,N-popSize)];
        end
    else
        p1 = p1(:)';
        if numel(p1) ~= N
            % pad or truncate to length N (non-interactive safe behavior)
            if numel(p1) > N
                p1 = p1(1:N);
            else
                p1 = [p1, randi(popSize,1,N-numel(p1))];
            end
        end
    end

    % ---- associate to reference directions (NSGA-III style) ----
    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end
    PopObj = Population.objs;
    PopObj = PopObj - repmat(Zmin,popSize,1);
    [~,M]  = size(PopObj);
    NZ     = size(Z,1);

    % Normalization (same as EnvironmentalSelection/LastSelection)
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for m = 1 : M
        [~,Extreme(m)] = min(max(PopObj./repmat(w(m,:),popSize,1),[],2));
    end
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a)) || any(isinf(a))
        a = max(PopObj,[],1)';
        a(a<=0) = 1;
    end
    PopObjN = PopObj./repmat(a',popSize,1);

    Cosine   = 1 - pdist2(PopObjN,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObjN.^2,2)),1,NZ).*sqrt(max(0,1-Cosine.^2));
    [~,pi]   = min(Distance',[],1);  % 1 x popSize

    state.pi     = pi(:);
    state.PopObj = PopObjN;

    % ---- select p2 for each p1 ----
    p2 = zeros(1,N);
    for t = 1:N
        i1 = p1(t);
        i1 = max(1,min(popSize,i1));
        niche1 = pi(i1);

        cand = find(pi ~= niche1);
        if isempty(cand)
            % only one niche present; fallback: random other individual
            if popSize == 1
                p2(t) = i1;
            else
                r = randi(popSize-1);
                if r >= i1, r = r + 1; end
                p2(t) = r;
            end
            continue;
        end

        switch mode
            case "random"
                p2(t) = cand(randi(numel(cand)));
            otherwise % "farthest"
                v1 = PopObjN(i1,:);
                V  = PopObjN(cand,:);
                d2 = sum((V - v1).^2,2);
                [~,idx] = max(d2);
                p2(t) = cand(idx);
        end
    end
end

