function intersects = check_line_obstacle_intersection(p1, p2, obstacle, obstacle_type)
% 统一的障碍物检测接口
% 输入:
%   p1, p2: 线段端点 [x, y]
%   obstacle: 障碍物数据
%   obstacle_type: 'rectangle' 或 'circle'
% 输出:
%   intersects: 布尔值，true表示相交，false表示不相交

    if strcmp(obstacle_type, 'rectangle')
        % 矩形：[x_min, y_min, x_max, y_max]
        intersects = check_line_rect_intersect(p1, p2, obstacle);
        
    elseif strcmp(obstacle_type, 'circle')
        % 圆形：[cx, cy, radius]
        intersects = check_line_circle_intersect(p1, p2, obstacle);
        
    else
        error('Unknown obstacle type: %s. Use ''rectangle'' or ''circle''.', obstacle_type);
    end
end

