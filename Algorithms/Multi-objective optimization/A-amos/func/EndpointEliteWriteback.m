function [Population,elite] = EndpointEliteWriteback(Population,elite)
%EndpointEliteWriteback - Maintain and write back endpoint elites.
%
%   [Population,elite] = EndpointEliteWriteback(Population,elite) maintains
%   two endpoint elites:
%     - elite.f1: solution with minimal objective 1
%     - elite.f2: solution with minimal objective 2
%   After updating elites from the current Population, if an elite is not
%   present in Population (by identical decision vector), it is written back
%   by replacing the worst solution in the corresponding objective (max f1
%   for elite.f1, max f2 for elite.f2).
%
%   This is a "zero extra FE" strategy: it reuses already-evaluated
%   solutions.

    if nargin < 2 || isempty(elite)
        elite = struct('f1',[],'f2',[]);
    end
    if isempty(Population) || isempty(Population.decs) || isempty(Population.objs)
        return;
    end

    Decs = Population.decs;
    Objs = Population.objs;
    M    = size(Objs,2);

    % Update elites from current population (lexicographic safety for M<2)
    [~,i1] = min(Objs(:,1));
    cand1  = Population(i1);
    elite  = updateElite(elite,'f1',cand1,1);
    if M >= 2
        [~,i2] = min(Objs(:,2));
        cand2  = Population(i2);
        elite  = updateElite(elite,'f2',cand2,2);
    end

    % Write back if missing
    Population = writebackIfMissing(Population,elite,'f1',1);
    if M >= 2
        Population = writebackIfMissing(Population,elite,'f2',2);
    end
end

function elite = updateElite(elite,field,cand,objIdx)
    if ~isfield(elite,field) || isempty(elite.(field))
        elite.(field) = cand;
        return;
    end
    bestObj = elite.(field).objs;
    candObj = cand.objs;
    if candObj(objIdx) < bestObj(objIdx)
        elite.(field) = cand;
    end
end

function Population = writebackIfMissing(Population,elite,field,objIdx)
    if ~isfield(elite,field) || isempty(elite.(field))
        return;
    end
    eliteSol = elite.(field);
    popDecs  = Population.decs;
    if isempty(popDecs)
        return;
    end
    present = any(all(popDecs == eliteSol.decs,2));
    if present
        return;
    end

    popObjs = Population.objs;
    if isempty(popObjs)
        return;
    end
    [~,worstIdx] = max(popObjs(:,objIdx));
    Population(worstIdx) = eliteSol;
end

