function [child1] = bianyi(child1)
    % 提取起点和终点的索引
    start_index_1 = find(child1 == 0);
    end_index_1 = find(child1 == -1);

   % 去掉起点和终点
    parent1_trimmed = [];
    for i = 1:length(start_index_1)
        parent1_trimmed = [parent1_trimmed, child1(start_index_1(i)+1:end_index_1(i)-1)];
    end
   
    crossover_point1 = randi([1, length(parent1_trimmed)-1]);
    crossover_point2 = randi([1, length(parent1_trimmed)-1]);
    
    % 确保两个交叉点不同
    while crossover_point1 == crossover_point2
        crossover_point2 = randi([1, length(parent1_trimmed)-1]);
    end
    
    d1=parent1_trimmed(crossover_point1);
    d2=parent1_trimmed(crossover_point2);
    parent1_trimmed(crossover_point1)=d2;
    parent1_trimmed(crossover_point2)=d1;

% 将起点和终点重新插入
    child1 = [];
    for i = 1:length(start_index_1)
        if i == 1
            child1 = [child1, 0, parent1_trimmed(1:end_index_1(i)-2), -1];
        else
        child1=[child1,0,parent1_trimmed(start_index_1(i)-2*(i-1):end_index_1(i)-2*(i-1)-1-1),-1];
        end
    end
   
   % 随机选择一段位置进行交叉
   random_index = randi([2, 3]);
   if start_index_1(random_index)-start_index_1(random_index-1)>3
      numindex=start_index_1(random_index)-2;
      num=child1(numindex);
      child1(numindex) = -1;
      child1(numindex+1) = 0;
      child1(numindex+2) = num;
   end
end