function [child1, child2] = genetic_crossover(parent1, parent2)
    % 提取起点和终点的索引
    start_index_1 = find(parent1 == 0);
    end_index_1 = find(parent1 == -1);
    start_index_2 = find(parent2 == 0);
    end_index_2 = find(parent2 == -1);
    
    % 去掉起点和终点
    parent1_trimmed = [];
    for i = 1:length(start_index_1)
        parent1_trimmed = [parent1_trimmed, parent1(start_index_1(i)+1:end_index_1(i)-1)];
    end
    parent2_trimmed = [];
    for i = 1:length(start_index_2)
        parent2_trimmed = [parent2_trimmed, parent2(start_index_2(i)+1:end_index_2(i)-1)];
    end
    
    % 随机选择一段位置进行交叉
    crossover_point = randi([1, length(parent1_trimmed)-1]);
    
    % 交叉操作
    child1_trimmed = [parent1_trimmed(1:crossover_point), parent2_trimmed(crossover_point+1:end)];
    child2_trimmed = [parent2_trimmed(1:crossover_point), parent1_trimmed(crossover_point+1:end)];
    
 % 修正重复基因
 % child1_trimmed=[2,3,4,4,2];
    child1_trimmed = fix_duplicates(child1_trimmed);
    child2_trimmed = fix_duplicates(child2_trimmed);



    % 将起点和终点重新插入
    child1 = [];
    
    for i = 1:length(start_index_1)
        if i == 1
            child1 = [child1, 0, child1_trimmed(1:end_index_1(i)-2), -1];
        else
        child1=[child1,0,child1_trimmed(start_index_1(i)-2*(i-1):end_index_1(i)-2*(i-1)-1-1),-1];
        end
    end
    child2 = [];
        for i = 1:length(start_index_2)
        if i == 1
            child2 = [child2, 0, child2_trimmed(1:end_index_2(i)-2), -1];
        else
            child2=[child2,0,child2_trimmed(start_index_2(i)-2*(i-1):end_index_2(i)-2*(i-1)-1-1),-1];
        end
    end
    
    % % 修正重复基因
    % child1 = fix_duplicates(child1, parent2);
    % child2 = fix_duplicates(child2, parent1);
end

% function child = fix_duplicates(child, parent)
%     % 找到重复基因的位置
%     [~, unique_idx] = unique(child);
%     duplicate_idx = setdiff(1:length(child), unique_idx);
%     % gene = 15-(sum(child)-child(i));   
%     % gene_idx_in_parent = find(parent == gene);
%     % % 替换重复基因
%     % child(i) = parent(gene_idx_in_parent);
% 
%     % 修正重复基因
%     for i = duplicate_idx
%         gene = child(i);
%         % 找到该基因在另一父代的对应位置
%         gene_idx_in_parent = find(parent == gene);
%         % 替换重复基因
%         child(i) = gene_idx_in_parent;
% 
%         % 如果继续出现重复值，则递归修正
%         child = fix_duplicates(child, parent);
%     end
% end
function child = fix_duplicates(child)
[B, unique_idx] = unique(child);
duplicate_idx = setdiff(1:length(child), unique_idx);
% 定义完整序列
full_sequence = 1:5;

% 找出缺失的元素
missing_elements = setdiff(full_sequence, B);
for i=1:size(duplicate_idx,2)
% 将缺失的元素添加到 B 中
child(duplicate_idx(i)) = missing_elements(i);
end
end