function Distance = CalFitness(PopObj,FrontNo, MaxFNo)

[N,M]        =  size(PopObj);
Distance     = zeros(1,N);

for i = 1 : MaxFNo
    Current      = find(FrontNo==i);
    PopC         = PopObj(Current,:);
    %% Calculate the hyperplane
    [~,Rank]     = sort(PopC,'descend');
    Extreme = FindCornerSolutions(PopC); 
    % Calculate the hyperplane
    Hyperplane = PopC(Extreme,:)\ones(length(Extreme),1);
    % Calculate the distance of each solution to the hyperplane
    Dist       = -(PopC*Hyperplane-1)./sqrt(sum(Hyperplane.^2));
    %% Calculate the projected SDE distance
    DistanceValue = F_distance(Dist,PopC);
    Distance(Current) = DistanceValue;
    clear DistanceValue
end

end