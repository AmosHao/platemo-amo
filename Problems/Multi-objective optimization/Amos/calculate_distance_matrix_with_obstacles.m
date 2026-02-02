function d_matrix = calculate_distance_matrix_with_obstacles(dotss, forbidden_zones, crowded_zones, noise_zones, penalty_forbidden, penalty_crowded, penalty_noise)
% 计算考虑障碍物的距离矩阵
% 输入:
%   dotss: 坐标矩阵 (n×3)，每行为 [x, y, z]
%   forbidden_zones: 矩形禁飞区矩阵 (m1×4)，每行为 [x_min, y_min, x_max, y_max]
%   crowded_zones: 圆形人流密集区矩阵 (m2×3)，每行为 [cx, cy, radius]
%   noise_zones: 矩形噪音敏感区矩阵 (m3×4)，每行为 [x_min, y_min, x_max, y_max]
%   penalty_forbidden: 禁飞区惩罚距离
%   penalty_crowded: 人流密集区惩罚距离
%   penalty_noise: 噪音敏感区惩罚距离
% 输出:
%   d_matrix: 距离矩阵 (n×n)

    n = size(dotss, 1);
    d_matrix = zeros(n, n);
    
    for i = 1:n
        for j = 1:n
            if i == j
                d_matrix(i, j) = 0;
            else
                p1 = dotss(i, 1:2);  % 起点坐标
                p2 = dotss(j, 1:2);  % 终点坐标
                
                % 计算欧氏距离
                d_euclidean = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
                
                % 检测是否经过障碍物，累加惩罚
                total_penalty = 0;
                
                % 检测禁飞区（矩形）
                if ~isempty(forbidden_zones)
                    for k = 1:size(forbidden_zones, 1)
                        if check_line_obstacle_intersection(p1, p2, ...
                            forbidden_zones(k, :), 'rectangle')
                            total_penalty = total_penalty + penalty_forbidden;
                            % 注意：如果一条线段经过多个障碍物，累加所有惩罚
                        end
                    end
                end
                
                % 检测人流密集区（圆形）
                if ~isempty(crowded_zones)
                    for k = 1:size(crowded_zones, 1)
                        if check_line_obstacle_intersection(p1, p2, ...
                            crowded_zones(k, :), 'circle')
                            total_penalty = total_penalty + penalty_crowded;
                        end
                    end
                end
                
                % 检测噪音敏感区（矩形）
                if ~isempty(noise_zones)
                    for k = 1:size(noise_zones, 1)
                        if check_line_obstacle_intersection(p1, p2, ...
                            noise_zones(k, :), 'rectangle')
                            total_penalty = total_penalty + penalty_noise;
                        end
                    end
                end
                
                % 最终距离 = 欧氏距离 + 惩罚
                d_matrix(i, j) = d_euclidean + total_penalty;
            end
        end
    end
end

