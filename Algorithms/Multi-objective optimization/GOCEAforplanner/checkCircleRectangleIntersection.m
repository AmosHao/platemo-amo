function isIntersect = checkCircleRectangleIntersection(r, circle_center, track_point1,track_point2, n2, rectangleVertices)
    % 计算两个法向量的夹角
    % n1 和 n2 是法向量
    % circleCenter 是圆心坐标
    % rectangleVertices 是长方形的四个顶点坐标
    % r 是圆的半径
    % 计算法向量
    vec1 = track_point1 - circle_center; % 圆心到第一航迹点的向量
    vec2 = track_point2 - circle_center; % 圆心到第二航迹点的向量
    n1 = cross(vec1, vec2);  % 使用叉乘得到法向量
    % 确保法向量是单位向量
    n1 = n1 / norm(n1);
    n2 = n2 / norm(n2);

    % 计算夹角
    cos_theta = dot(n1, n2);
    % 确保夹角计算在合法范围内
    cos_theta = min(1, max(-1, cos_theta));
    theta = acos(cos_theta);
    
    % 计算 d1 = r * sin(theta)
    d1 = r * sin(theta);

    % 计算圆心到长方形平面的垂直距离 d2
    % 假设长方形平面可以通过长方形的顶点和法向量得到
    % 使用顶点的第一个点作为平面上的点 (x0, y0, z0)
    x0 = rectangleVertices(1, 1);
    y0 = rectangleVertices(1, 2);
    z0 = rectangleVertices(1, 3);

    % 平面方程 Ax + By + Cz + D = 0
    % 通过法向量 n2 确定平面方程
    A = n2(1);
    B = n2(2);
    C = n2(3);
    D = -(A * x0 + B * y0 + C * z0);

    % 圆心到平面的垂直距离 d2
    circleCenter = circleCenter(:)'; % 确保是行向量
    d2 = abs(A * circleCenter(1) + B * circleCenter(2) + C * circleCenter(3) + D) / norm([A, B, C]);

    % 判断是否相交
    isIntersect = d1 > d2;
end

