function intersects = check_line_rect_intersect(p1, p2, rect)
% 检测线段是否与轴对齐矩形相交
% 输入:
%   p1, p2: 线段两个端点 [x, y]
%   rect: 矩形 [x_min, y_min, x_max, y_max]
% 输出:
%   intersects: 布尔值，true表示相交，false表示不相交

    % 方法1：快速拒绝测试（Fast Rejection Test）
    % 如果线段完全在矩形的一侧，则不相交
    x_min = min(p1(1), p2(1));
    x_max = max(p1(1), p2(1));
    y_min = min(p1(2), p2(2));
    y_max = max(p1(2), p2(2));
    
    % 如果线段的包围盒与矩形不相交，则线段不相交
    if x_max < rect(1) || x_min > rect(3) || ...
       y_max < rect(2) || y_min > rect(4)
        intersects = false;
        return;
    end
    
    % 方法2：检查线段端点是否在矩形内
    if point_in_rect(p1, rect) || point_in_rect(p2, rect)
        intersects = true;
        return;
    end
    
    % 方法3：检测线段是否与矩形的任意一条边相交
    % 使用参数方程：p = p1 + t*(p2-p1), t in [0,1]
    
    % 计算线段方向向量
    dx = p2(1) - p1(1);
    dy = p2(2) - p1(2);
    
    % 检测与矩形四条边的交点
    % 左边界：x = rect(1)
    if abs(dx) > 1e-10  % 避免除零
        t = (rect(1) - p1(1)) / dx;
        if t >= 0 && t <= 1
            y = p1(2) + t * dy;
            if y >= rect(2) && y <= rect(4)
                intersects = true;
                return;
            end
        end
    end
    
    % 右边界：x = rect(3)
    if abs(dx) > 1e-10
        t = (rect(3) - p1(1)) / dx;
        if t >= 0 && t <= 1
            y = p1(2) + t * dy;
            if y >= rect(2) && y <= rect(4)
                intersects = true;
                return;
            end
        end
    end
    
    % 下边界：y = rect(2)
    if abs(dy) > 1e-10
        t = (rect(2) - p1(2)) / dy;
        if t >= 0 && t <= 1
            x = p1(1) + t * dx;
            if x >= rect(1) && x <= rect(3)
                intersects = true;
                return;
            end
        end
    end
    
    % 上边界：y = rect(4)
    if abs(dy) > 1e-10
        t = (rect(4) - p1(2)) / dy;
        if t >= 0 && t <= 1
            x = p1(1) + t * dx;
            if x >= rect(1) && x <= rect(3)
                intersects = true;
                return;
            end
        end
    end
    
    intersects = false;
end

function inside = point_in_rect(p, rect)
% 检查点是否在矩形内
    inside = (p(1) >= rect(1) && p(1) <= rect(3) && ...
              p(2) >= rect(2) && p(2) <= rect(4));
end

