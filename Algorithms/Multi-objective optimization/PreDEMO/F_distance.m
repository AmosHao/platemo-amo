function DistanceValue = F_distance(Dist,FunctionValue)

    [N,M] = size(FunctionValue);
    PopObj = FunctionValue;
%%  sumof OBJECTIVE
     [~,rank] = sort(Dist,'descend');
    
 %%%%%%%%%%%%% % SDE with Sum of Objectives  %%%%%%%%%%%%%%%%%%%%
    DistanceValue = zeros(1,N);   
  
    for j = 2 : N
        SFunctionValue = max(PopObj(rank(1:j-1),:),repmat(PopObj(rank(j),:),(j-1),1));
        Distance = inf(1,j-1);   
        for i = 1 : (j-1)
            Distance(i) = norm(SFunctionValue(i,:)-PopObj(rank(j),:))/M;
        end           
        Distance               = min(Distance); 
        DistanceValue(rank(j)) = exp(-Distance);  
    end  
end