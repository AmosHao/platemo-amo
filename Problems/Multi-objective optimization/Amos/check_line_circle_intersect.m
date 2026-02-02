function intersects = check_line_circle_intersect(p1, p2, circle)
% 检测线段是否与圆形相交
% 输入:
%   p1, p2: 线段两个端点 [x, y]
%   circle: 圆形 [cx, cy, radius]
% 输出:
%   intersects: 布尔值，true表示相交，false表示不相交

    cx = circle(1);
    cy = circle(2);
    r = circle(3);
    
    % 方法：计算圆心到线段的最短距离
    % 如果最短距离 <= 半径，则相交
    
    % 计算线段方向向量
    dx = p2(1) - p1(1);
    dy = p2(2) - p1(2);
    len_sq = dx^2 + dy^2;
    
    if len_sq < 1e-10
        % 线段退化为点
        dist = sqrt((p1(1) - cx)^2 + (p1(2) - cy)^2);
        intersects = (dist <= r);
        return;
    end
    
    % 计算从p1到圆心的向量
    to_circle_x = cx - p1(1);
    to_circle_y = cy - p1(2);
    
    % 计算投影参数 t（在线段上的位置）
    % t = dot((circle - p1), (p2 - p1)) / ||p2 - p1||^2
    t = (to_circle_x * dx + to_circle_y * dy) / len_sq;
    
    % 将 t 限制在 [0, 1] 范围内
    t = max(0, min(1, t));
    
    % 计算线段上距离圆心最近的点
    closest_x = p1(1) + t * dx;
    closest_y = p1(2) + t * dy;
    
    % 计算最近点到圆心的距离
    dist = sqrt((closest_x - cx)^2 + (closest_y - cy)^2);
    
    intersects = (dist <= r);
end

