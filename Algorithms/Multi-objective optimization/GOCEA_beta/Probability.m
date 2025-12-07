function [BETA, L, K] = Probability(pop_b, pop, H, L, Kmax, t, Tu) 
    si = ones(size(pop, 1), 1); % 存活标志 si 
    % for i = 1:size(pop, 1)
    %     if ~isequal(pop(i, :), pop_b(i, :))
    %         si(i) = 0;
    %     end
    % end
    for i = 1:size(pop, 1)
        for j = 1:size(pop_b, 1)
            if ismember(pop(i, :), pop_b, 'rows')
                si(i) = 0; % 如果当前个体包含于上一代中，则将存活标志设置为0
                break; 
            end
        end
    end

    for i = 1:size(pop, 1) % 计算存活长度 L
        if si(i) == 0
            L(i) = min(L(i) + 1, H);
        else
            L(i) = 1;
        end
    end

    % La=max(L(i));
    % Lb=min(L(i));

    La=max(L(:,1));
    Lb=min(L(:,1));
    for i = 1:size(pop, 1)
       BETA(i) = (L(i) - Lb) / (La - Lb);% 更新交配限制概率 beta
    end
    
    if mod(t, Tu) == 0   % 更新最大聚类数 K
        K = ceil(Kmax*mean(BETA,1));
    else
        K = Kmax; % 不需要更新 K
    end
end
