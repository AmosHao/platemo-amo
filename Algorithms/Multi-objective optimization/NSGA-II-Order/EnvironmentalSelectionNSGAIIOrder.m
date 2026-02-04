function [Population,FrontNo,CrowdDis] = EnvironmentalSelectionNSGAIIOrder(Population,N)
% Environmental selection of NSGA-II with enhanced diversity for identical objectives
% 增强版环境选择：处理f1值相同的情况

    %% Non-dominated sorting
    PopObj = Population.objs;
    [FrontNo,MaxFNo] = NDSort(PopObj,length(Population));
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(PopObj,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% 增强：处理f1值相同的情况
    % 如果第一前沿中有很多f1值相同的解，使用决策空间距离作为辅助指标
    Front1Indices = find(FrontNo == 1);
    if length(Front1Indices) > N * 0.5  % 如果第一前沿解太多
        Front1Objs = PopObj(Front1Indices, :);
        Front1Decs = Population(Front1Indices).decs;
        
        % 检查f1值的唯一性
        uniqueF1 = unique(Front1Objs(:, 1));
        if length(uniqueF1) < length(Front1Indices) * 0.5  % 如果f1值多样性不足
            % 计算决策空间距离矩阵
            decDistances = pdist2(Front1Decs, Front1Decs);
            decDistances(logical(eye(size(decDistances)))) = Inf;
            
            % 为每个解计算最小决策空间距离（多样性指标）
            minDecDist = min(decDistances, [], 2);
            
            % 结合目标空间和决策空间的多样性
            combinedDiversity = CrowdDis(Front1Indices) + 0.1 * minDecDist / (max(minDecDist) + 1e-10);
            
            % 重新排序
            [~, newRank] = sort(combinedDiversity, 'descend');
            Front1Indices = Front1Indices(newRank);
            
            % 更新Next
            Next = false(size(Next));
            Next(Front1Indices(1:min(N, length(Front1Indices)))) = true;
            
            % 如果还需要更多解，从其他前沿选择
            if sum(Next) < N
                remaining = find(~Next);
                [~, rank] = sort(CrowdDis(remaining), 'descend');
                Next(remaining(rank(1:N-sum(Next)))) = true;
            end
        end
    end
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

function CrowdDis = CrowdingDistance(PopObj,FrontNo)
% Calculate the crowding distance of each solution

    [N,M] = size(PopObj);
    CrowdDis = zeros(1,N);
    Fronts = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontNo==Fronts(f));
        Fmax = max(PopObj(Front,:),[],1);
        Fmin = min(PopObj(Front,:),[],1);
        for i = 1 : M
            [~,Rank] = sortrows(PopObj(Front,i));
            CrowdDis(Front(Rank(1)))   = inf;
            CrowdDis(Front(Rank(end))) = inf;
            for j = 2 : length(Front)-1
                CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
            end
        end
    end
end
