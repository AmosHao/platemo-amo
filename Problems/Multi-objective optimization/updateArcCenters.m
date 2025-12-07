function PopDec = updateArcCenters(PopDec)
    % PopDec: 原始矩阵
    % r: 圆弧的半径
    
    [numIndividuals, numVars] = size(PopDec);

    n=(numVars-3)/6;
    % if numVars ~= 63
    %     error('PopDec的列数必须为63。');
    % end

    % 起点和终点
    startPoint = [0, 0, 0];
    endPoint = [2, 2, 0.3];

    % 处理每一行
    for i = 1:numIndividuals
        % 提取航迹点坐标
        x = PopDec(i, 1:n); % 航迹点x坐标
        y = PopDec(i, n+1:2*n); % 航迹点y坐标
        z = PopDec(i, 2*n+1:3*n); % 航迹点z坐标

        % 初始化圆心坐标矩阵
        arcCenters = zeros(n+1, 3); % 11个圆心坐标，每个坐标3个值

        % 计算圆心坐标
        for j = 0:10
            if j == 0
                % 第一个航迹点与起点生成圆心坐标
                p1 = startPoint;
                p2 = [x(j+1), y(j+1), z(j+1)];
            elseif j == 10
                % 第十个航迹点与终点生成圆心坐标
                p1 = [x(j), y(j), z(j)];
                p2 = endPoint;
            else
                % 中间的航迹点与相邻点生成圆心坐标
                p1 = [x(j), y(j), z(j)];
                p2 = [x(j+1), y(j+1), z(j+1)];
            end
            
            % 计算中点
            midpoint = (p1 + p2) / 2;
            
            % 计算两个点之间的距离
            distance = norm(p1 - p2);
            r=distance/2+0.05;
            % 计算圆心的偏移量
            if distance < 2 * r
                % 假设圆心在与当前点和目标点垂直的方向上
                offset = sqrt(r^2 - (distance / 2)^2); % 圆心到中点的距离
                direction = cross(p2 - p1, [0, 0, 1]); % 计算方向向量
                
                if norm(direction) ~= 0
                    direction = direction / norm(direction); % 单位化方向向量
                    arc_center = midpoint + offset * direction;
                    
                    % 存储圆心坐标
                    arcCenters(j+1, :) = arc_center;
                end
            end
        end

        % 存储结果到 PopDec
PopDec(i, 3*n+1:4*n+1) = arcCenters(:, 1)';  % 第一列写入 31-41
PopDec(i, 4*n+2:5*n+2) = arcCenters(:, 2)';  % 第二列写入 42-52
PopDec(i, 5*n+3:6*n+3) = arcCenters(:, 3)';  % 第三列写入 53-63
    end
end