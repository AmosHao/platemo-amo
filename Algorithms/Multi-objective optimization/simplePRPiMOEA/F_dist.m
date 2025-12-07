function D = F_dist(PopObj)

%% Calculate D(i)
   [N,~]  = size(PopObj);
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D = 1./(Distance(:,floor(sqrt(N)))+2);
    
end