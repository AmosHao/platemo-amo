% pop=[3,3,1.5,4,4,1.5,3.5,3.5,1.5];
% cubesdata=[1,1,0,2,1,0,2,2,0,1,2,0,1,1,3,2,1,3,2,2,0,1,2,3];
% varDim=21;
% col=check_coll(pop, cubesdata,varDim)

% filename = 'pop2.xlsx';
% sheet = 1;  % 假设数据在第一个工作表中
% % 使用xlsread读取数据
% pop = xlsread(filename, sheet);
% % 指定Excel文件名和工作表
% filename = 'cubes_data_curve_small_2.xlsx';
% sheet = 1;  % 假设数据在第一个工作表中
% % 使用xlsread读取数据
% cubesdata = xlsread(filename, sheet);
% varDim
% collision_indices = check_collision(pop, cubesdata,varDim);

function collision_indices = check_collisions(pop, cubesdata,varDim)
    % 输入参数：
    % pop - 200x63的矩阵，表示200条航迹线的数据
    % cubesdata - 5x24的矩阵，表示5个长方体的八个顶点坐标
    % writematrix(pop, 'pop2.xlsx', 'Sheet', 'dec_data');
    % 定义长方体数量和航迹线数量
    num_cubes = size(cubesdata, 1);
    num_trajectories = size(pop, 1);

    % 初始化碰撞结果的索引
    collision_indices = false(num_trajectories, 1);
    n=(varDim-3)/6;
    % 遍历每条航迹线
    for t = 1:num_trajectories
        % 获取当前航迹线的圆数据
        x = pop(t, 1:n);
        y = pop(t, n+1:2*n);
        z = pop(t, 2*n+1:3*n);
        xStart=[0,0,0]; xEnd=[2,2,0.3];                              % 起点和终点
        x=[xStart(1),x,xEnd(1)]; y=[xStart(2),y,xEnd(2)]; z=[xStart(3),z,xEnd(3)];
        center_x = pop(t, 3*n+1:4*n+1);
        center_y = pop(t, 4*n+2:5*n+2);
        center_z = pop(t, 5*n+3:6*n+3);

        
        % 计算每个圆的半径
        radii = calculate_radii(x, y, z,center_x,center_y,center_z); 
        collision_yuan = false(n+1, num_cubes);

        for j=1:n+1%遍历11个圆
            circle_center = [center_x(j), center_y(j), center_z(j)];
            
        % 遍历每个长方体
        for c = 1:num_cubes
            % 获取当前长方体的8个顶点
            cube_vertices = reshape(cubesdata(c, :), [3, 8]);
            
            % 计算长方体的6个面
            faces= [
                cube_vertices(:,[3, 2, 1, 4]); % 第一个面的四个顶点
                cube_vertices(:,[5, 6, 7, 8]);  % 第二个面的四个顶点
                cube_vertices(:,[5, 1, 2, 6]);  % 第三个面的四个顶点
                cube_vertices(:,[6, 2, 3, 7]);  % 第四个面的四个顶点
                cube_vertices(:,[7, 3, 4, 8]);  % 第五个面的四个顶点
                cube_vertices(:,[8, 4, 1, 5]);  % 第六个面的四个顶点
            ];
            %计算长方体每个面的法向量和中心点
            face_normals = zeros(3, 6);
            face_centers = zeros(3, 6);
            for i = 1:6
                ii = (i-1) * 3 + 1;
                face_vertices = faces(ii:ii+2,:);
                [normal, center] = compute_face_normal_and_center(face_vertices);%输入每个面四个点坐标
                face_normals(:, i) = normal;
                face_centers(:, i) = center;
            end

 
                % 检查圆是否与长方体的任何一个面相交
                collision_found = true(1,6);%6个面都是true则collision_yuan（j，c）是true
                for k = 1:6
                kk = (k-1) * 3 + 1;
                    face_vertices = faces(kk:kk+2,:);
                    if check_circle_intersects_face(circle_center, face_normals(:, k), face_centers(:, k), face_vertices,radii(j))%点投影在面内，用method进一步判断
                        track_point1=[x(j), y(j), z(j)];
                        track_point2=[x(j+1), y(j+1), z(j+1)];
                        if checkCircleRectangleIntersection(radii(j), circle_center, track_point1,track_point2, face_normals(:, k), face_vertices)
                        collision_found(1,k) = false;
                        end
                    end
                end

                if all(collision_found)
                collision_yuan(j, c) = true;
                end      
                

            %检查圆是否包围长方体
            
            radius = radii(j);
            track_point1=[x(j), y(j), z(j)];
            track_point2=[x(j+1), y(j+1), z(j+1)];
            isContained = check_circle_contains_box(circle_center, radius, track_point1, track_point2, cube_vertices');
            if isContained
            collision_yuan(j,c) = isContained;
            end
                

          end
          end
    if all(collision_yuan)
        collision_indices(t) = true;
    end
    end
end

% 计算每个圆的半径
function radii = calculate_radii(x, y, z,center_x,center_y,center_z)
    num_points = size(center_x, 2);
    radii = zeros(1, num_points);
    
    for i = 1:num_points
        radii(i) = sqrt((x(i)-center_x(i)).^2+(y(i)-center_y(i)).^2+(z(i)-center_z(i)).^2);
    end
end

% 计算长方体面法向量和中心点
function [normal, center] = compute_face_normal_and_center(vertices)
    p1 = vertices(:,1);
    p2 = vertices(:,2);
    p3 = vertices(:,3);
    
    v1 = p2 - p1;
    v2 = p3 - p2;
    normal = cross(v1, v2);
    normal = normal / norm(normal);
    
    center = (p1 + p3) / 2;
end


% function is_inside = check_circle_intersects_face(circle_center, face_normal, face_center, face_vertices,circle_radius)
%     % 计算面与圆心的投影点
%     d = dot(face_normal, (circle_center' - face_center));
%     proj_point = circle_center' - d * face_normal;
% 
%     % 投影点坐标
%     px = proj_point(1);
%     py = proj_point(2);
%     pz = proj_point(3);
% 
%     % 矩形顶点坐标
%     face_vertices = face_vertices';
%     x_coords = face_vertices(:, 1);
%     y_coords = face_vertices(:, 2);
%     z_coords = face_vertices(:, 3);
% 
%     % % 计算矩形对角线的长度
%     % diagonal_length = sqrt((max(x_coords) - min(x_coords))^2 + ...
%     %                        (max(y_coords) - min(y_coords))^2 + ...
%     %                        (max(z_coords) - min(z_coords))^2);
% 
%     % 计算投影点到矩形四个顶点的距离
%     distances = sqrt((px - x_coords).^2 + ...
%                      (py - y_coords).^2 + ...
%                      (pz - z_coords).^2);
% 
% 
%     % 计算最小和最大坐标
%     min_x = min(x_coords);
%     max_x = max(x_coords);
%     min_y = min(y_coords);
%     max_y = max(y_coords);
% 
%     % 检查投影点是否在矩形内部
%     if min_x <= px && px <= max_x && min_y <= py && py <= max_y
%         is_inside = (min_x <= px && px <= max_x && min_y <= py && py <= max_y);
%     % 检查是否距离最大值小于等于对角线长度加上半径
%     else 
%         is_inside = min(distances) >= circle_radius;
%     end
% end
function is_inside = check_circle_intersects_face(circle_center, face_normal, face_center, face_vertices, circle_radius)
    % 计算面与圆心的投影点
    d = dot(face_normal, (circle_center' - face_center));
    proj_point = circle_center' - d * face_normal;
    
    % 投影点坐标
    px = proj_point(1);
    py = proj_point(2);
    
    % 矩形顶点坐标
    face_vertices = face_vertices';
    x_coords = face_vertices(:, 1);
    y_coords = face_vertices(:, 2);
    
    % 计算矩形的边界
    min_x = min(x_coords);
    max_x = max(x_coords);
    min_y = min(y_coords);
    max_y = max(y_coords);

    % 检查投影点是否在矩形内部
    if min_x <= px && px <= max_x && min_y <= py && py <= max_y
is_inside = (min_x <= px && px <= max_x && min_y <= py && py <= max_y);
    else
        % 投影点在矩形外部，检查圆心到矩形四个顶点的距离
        distances = sqrt((px - x_coords).^2 + (py - y_coords).^2);
        is_inside = min(distances) <= circle_radius*0.0001;
    end
end
% % 判断投影点是否在多边形内
% function is_inside = point_in_polygon(point, polygon_vertices)
%     % 使用射线法判断点是否在多边形内
%     n = size(polygon_vertices, 2); % 多边形顶点数
%     is_inside = false;
% 
%     % 遍历每条边
%     j = n;
%     for i = 1:n
%         xi = polygon_vertices(1, i);
%         yi = polygon_vertices(2, i);
%         xj = polygon_vertices(1, j);
%         yj = polygon_vertices(2, j);
% 
%         % 判断点是否在多边形边的上方
%         if ((yi > point(2)) ~= (yj > point(2))) && ...
%            (point(1) < (xj - xi) * (point(2) - yi) / (yj - yi) + xi)
%             is_inside = ~is_inside;
%         end
% 
%         j = i;
%     end
% end

function isContained = check_circle_contains_box(circle_center, radius, track_point1, track_point2, box_vertices)
    % circle_center: 圆心坐标 (1x3)
    % radius: 圆的半径
    % track_point1: 第一航迹点 (1x3)
    % track_point2: 第二航迹点 (1x3)
    % box_vertices: 长方体的顶点坐标 (Nx3)

    % 计算法向量
    vec1 = track_point1 - circle_center; % 圆心到第一航迹点的向量
    vec2 = circle_center-track_point2; % 圆心到第二航迹点的向量
    normal_vector = cross(vec1, vec2);  % 使用叉乘得到法向量

    % 计算长方体的投影点到圆所在平面
    projected_points = zeros(size(box_vertices)); % 初始化投影点矩阵
    for i = 1:size(box_vertices, 1)
        % 计算当前顶点到圆心的向量
        point = box_vertices(i, :);
        vector_to_circle_center = point - circle_center;
        
        % 计算点到平面的投影
        distance_to_plane = dot(vector_to_circle_center, normal_vector) / norm(normal_vector);
        projected_point = point - (distance_to_plane * normal_vector / norm(normal_vector));
        
        projected_points(i, :) = projected_point; % 存储投影点
    end

    % 计算投影点到圆心的距离
    distances = sqrt(sum((projected_points - circle_center).^2, 2));

    % 检查是否所有投影点的距离都大于圆的半径
    if all(distances < radius)
        isContained = true; % 圆包围长方体
    else
        isContained = false; % 圆与长方体可能相交
    end
end

function isIntersect = checkCircleRectangleIntersection(r, circle_center, track_point1,track_point2, n2, rectangleVertices)
    % 计算两个法向量的夹角
    % n1 和 n2 是法向量
    % circleCenter 是圆心坐标
    % rectangleVertices 是长方形的四个顶点坐标
    % r 是圆的半径
    % 计算法向量
    vec1 = track_point1 - circle_center; % 圆心到第一航迹点的向量
    vec2 = circle_center-track_point2; % 圆心到第二航迹点的向量
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
    circle_center = circle_center(:)'; % 确保是行向量
    d2 = abs(A * circle_center(1) + B * circle_center(2) + C * circle_center(3) + D) / norm([A, B, C]);

    % 判断是否相交
    isIntersect = d1 > d2;
end

