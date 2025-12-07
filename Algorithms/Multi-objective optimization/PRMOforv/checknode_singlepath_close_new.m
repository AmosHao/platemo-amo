function [currentSol,validPoints]=checknode_singlepath_close_new(currentSol,varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2)
%读取cubes数据
 % 文件名
% filename = 'city_environment_data_6.xlsx';
% % 读取长方体数据
% buildings = xlsread(filename, 'Buildings');
% numBuildings = size(buildings, 1);
% % 读取圆柱体数据
% cylinders = xlsread(filename, 'Cylinders');
% numCylinders = size(cylinders, 1);
% % 读取四棱锥数据
% pyramids = xlsread(filename, 'Pyramids');
% numPyramids = size(pyramids, 1);
% % 读取球体数据
% spheres = xlsread(filename, 'Spheres');
% numSpheres = size(spheres, 1);

% pop=xlsread('travel02', 'dec_data');
n_dots=varDim/4;
% x=[min(dot1(1),dot2(1)),currentSol(:,1:n_dots),max(dot1(1),dot2(1))];
% y=[min(dot1(2),dot2(2)),currentSol(:,n_dots+1:2*n_dots),max(dot1(2),dot2(2))];
% z=[min(dot1(3),dot2(3)),currentSol(:,2*n_dots+1:3*n_dots),max(dot1(3),dot2(3))];
x=[dot1(1),currentSol(:,1:n_dots),dot2(1)];
y=[dot1(2),currentSol(:,n_dots+1:2*n_dots),dot2(2)];
z=[dot1(3),currentSol(:,2*n_dots+1:3*n_dots),dot2(3)];


validPoints = false(n_dots, 1);
% 设置距离阈值（可以根据需要进行调整）
distanceThreshold = 1400;  % 例如：5单位

for i = 1:n_dots
    % if i==7
    % aaa=113
    % end
    % 当前候选点坐标
    candidateCoords = ([x(:, i), y(:, i), z(:, i)]+[x(:, i+2), y(:, i+2), z(:, i+2)])/2;
    % 检查该点是否与任何一个长方体发生碰撞
    isValid = true; % 默认认为是有效（在外部）

    % 检查球体
    % tic;  % 开始计时
    for k = 1:numSpheres
        % 获取球体参数
        center = spheres(k, 1:3);
        radius = spheres(k, 4);
        
        % 计算点到球体圆心的距离
        distToCenter = norm(candidateCoords - center);
        
        % 只有在距离较近的情况下检查碰撞
        % if distToCenter <= (radius + distanceThreshold)
            % 检查点是否在球体内
            if distToCenter <= radius
                isValid = false;
                break; % 退出当前球体的检测
            end
        % end
    end
    % toc;


    if isValid
        % 检查圆柱体
        % tic;
        for k = 1:numCylinders
            % 获取圆柱体参数
            % if k==117
            %     bb=1111;
            % end
            center = cylinders(k, 1:2);
            radius = cylinders(k, 3);
            height = cylinders(k, 4);
            
            % 计算点到圆柱体圆心的距离
            distToCenter = norm(candidateCoords(1:2) - center(1:2));
            
            % 只有在距离较近的情况下检查碰撞
            % if distToCenter <= (radius + distanceThreshold) && candidateCoords(3) <= height
                if distToCenter <= radius && candidateCoords(3) <= height
                    isValid = false;
                    break; % 退出当前圆柱体的检测
                end
            % end
        end
    end
    % toc;
    if isValid
        % 检查四棱锥
        % tic;
        for k = 1:numPyramids
            % 获取四棱锥的5个顶点
            vertices_p = pyramids(k, :);
            vertices_p = reshape(vertices_p, 3, 5);
            

            
            % 计算每个顶点的距离
            distToPyramid = min(vecnorm(vertices_p - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToPyramid <= distanceThreshold
                            % 计算四棱锥的面
            faces_p = [
                vertices_p(:,[3, 2, 1]); % 底面
                vertices_p(:,[1,5,4]);  % 第二个面的四个顶点
                vertices_p(:,[4,5,3]);  % 第三个面的四个顶点
                vertices_p(:,[2,5,1]);  % 第四个面的四个顶点
                vertices_p(:,[3,5,2]);  % 第五个面的四个顶点
            ];
                countInside_p = 0;
                for j = 1:3:15
                    % 面的三个顶点
                    v1_p = faces_p(j:j+2, 1);
                    v2_p = faces_p(j:j+2, 2);
                    v3_p = faces_p(j:j+2, 3);
                    % 计算面的法向量
                    normal_p = cross(v2_p - v1_p, v3_p - v2_p);
                    normal_p = normal_p / norm(normal_p); % 单位法向量
                    % 计算点到面的向量
                    vector_p = candidateCoords - v1_p';
                    % 计算点到面的向量与法向量的点积
                    dotProduct_p = dot(vector_p, normal_p');
                    % 如果点在面的内部，则增加内部计数
                    if dotProduct_p <= 0
                        countInside_p = countInside_p + 1;
                    end
                end
                % 如果点在所有面的内部，标记为无效
                if countInside_p == 5
                    isValid = false;
                    break; % 退出当前四棱锥的检测
                end
            end
        end
    end
    % toc;

     if isValid
        % 检查长方体
        % tic;
        % for k = 1:numjinfei
            % 获取长方体的8个顶点
            vertices = jinfei(1, :);
            % 将 vertices 重新排列成 8x3 的矩阵
            vertices_matrix = reshape(vertices, 3, 8);
           
            
            % 计算长方体的距离
            distToBuilding = min(vecnorm(vertices_matrix - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToBuilding <= distanceThreshold
                 % 定义长方体的面
            faces = [
                vertices_matrix(:,[3, 2, 1, 4]); % 第一个面的四个顶点
                vertices_matrix(:,[5, 6, 7, 8]);  % 第二个面的四个顶点
                vertices_matrix(:,[5, 1, 2, 6]);  % 第三个面的四个顶点
                vertices_matrix(:,[6, 2, 3, 7]);  % 第四个面的四个顶点
                vertices_matrix(:,[7, 3, 4, 8]);  % 第五个面的四个顶点
                vertices_matrix(:,[8, 4, 1, 5]);  % 第六个面的四个顶点
            ];
                countInside = 0;
                for j = 1:3:18
                    % 面的三个顶点
                    v1 = faces(j:j+2, 1);
                    v2 = faces(j:j+2, 2);
                    v3 = faces(j:j+2, 3);
                    % 计算面的法向量
                    normal = cross(v2 - v1, v3 - v2);
                    normal = normal / norm(normal); % 单位法向量
                    % 计算点到面的向量
                    vector = candidateCoords - v1';
                    % 计算点到面的向量与法向量的点积
                    dotProduct = dot(vector, normal');
                    % 如果点在面的内部，则增加内部计数
                    if dotProduct <= 0
                        countInside = countInside + 1;
                    end
                end
                % 如果点在所有六个面的内部，标记为无效
                if countInside == 6
                    isValid = false;
                    % break; % 退出当前长方体的检测
                end
            end
        % end
     end
     % toc;
    if isValid
        % 检查长方体
        % tic;
        for k = 1:numBuildings
            % 获取长方体的8个顶点
            vertices = buildings(k, :);
            % 将 vertices 重新排列成 8x3 的矩阵
            vertices_matrix = reshape(vertices, 3, 8);
           
            
            % 计算长方体的距离
            distToBuilding = min(vecnorm(vertices_matrix - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToBuilding <= distanceThreshold
                 % 定义长方体的面
            faces = [
                vertices_matrix(:,[3, 2, 1, 4]); % 第一个面的四个顶点
                vertices_matrix(:,[5, 6, 7, 8]);  % 第二个面的四个顶点
                vertices_matrix(:,[5, 1, 2, 6]);  % 第三个面的四个顶点
                vertices_matrix(:,[6, 2, 3, 7]);  % 第四个面的四个顶点
                vertices_matrix(:,[7, 3, 4, 8]);  % 第五个面的四个顶点
                vertices_matrix(:,[8, 4, 1, 5]);  % 第六个面的四个顶点
            ];
                countInside = 0;
                for j = 1:3:18
                    % 面的三个顶点
                    v1 = faces(j:j+2, 1);
                    v2 = faces(j:j+2, 2);
                    v3 = faces(j:j+2, 3);
                    % 计算面的法向量
                    normal = cross(v2 - v1, v3 - v2);
                    normal = normal / norm(normal); % 单位法向量
                    % 计算点到面的向量
                    vector = candidateCoords - v1';
                    % 计算点到面的向量与法向量的点积
                    dotProduct = dot(vector, normal');
                    % 如果点在面的内部，则增加内部计数
                    if dotProduct <= 0
                        countInside = countInside + 1;
                    end
                end
                % 如果点在所有六个面的内部，标记为无效
                if countInside == 6
                    isValid = false;
                    break; % 退出当前长方体的检测
                end
            end
        end
    end
    % toc;

    % 如果点在所有形状的外部，标记为有效
    if isValid
        validPoints(i) = true;
        currentSol(:,i)=candidateCoords(1,1);
        currentSol(:,i+n_dots)=candidateCoords(1,2);
        currentSol(:,i+2*n_dots)=candidateCoords(1,3);
    else%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        candidateCoords = [x(:, i+1), y(:, i+1), z(:, i+1)];
        isValid = true; % 默认认为是有效（在外部）
        for k = 1:numSpheres
        % 获取球体参数
        center = spheres(k, 1:3);
        radius = spheres(k, 4);
        
        % 计算点到球体圆心的距离
        distToCenter = norm(candidateCoords - center);
        
        % 只有在距离较近的情况下检查碰撞
        % if distToCenter <= (radius + distanceThreshold)
            % 检查点是否在球体内
            if distToCenter <= radius
                isValid = false;
                break; % 退出当前球体的检测
            end
        % end
    end
    % toc;


    if isValid
        % 检查圆柱体
        % tic;
        for k = 1:numCylinders
            % 获取圆柱体参数
            % if k==117
            %     bb=11112;
            % end
            center = cylinders(k, 1:2);
            radius = cylinders(k, 3);
            height = cylinders(k, 4);
            
            % 计算点到圆柱体圆心的距离
            distToCenter = norm(candidateCoords(1:2) - center(1:2));
            
            % 只有在距离较近的情况下检查碰撞
            % if distToCenter <= (radius + distanceThreshold) && candidateCoords(3) <= height
                if distToCenter <= radius && candidateCoords(3) <= height
                    isValid = false;
                    break; % 退出当前圆柱体的检测
                end
            % end
        end
    end
    % toc;
    if isValid
        % 检查四棱锥
        % tic;
        for k = 1:numPyramids
            % 获取四棱锥的5个顶点
            vertices_p = pyramids(k, :);
            vertices_p = reshape(vertices_p, 3, 5);
            

            
            % 计算每个顶点的距离
            distToPyramid = min(vecnorm(vertices_p - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToPyramid <= distanceThreshold
                            % 计算四棱锥的面
            faces_p = [
                vertices_p(:,[3, 2, 1]); % 底面
                vertices_p(:,[1,5,4]);  % 第二个面的四个顶点
                vertices_p(:,[4,5,3]);  % 第三个面的四个顶点
                vertices_p(:,[2,5,1]);  % 第四个面的四个顶点
                vertices_p(:,[3,5,2]);  % 第五个面的四个顶点
            ];
                countInside_p = 0;
                for j = 1:3:15
                    % 面的三个顶点
                    v1_p = faces_p(j:j+2, 1);
                    v2_p = faces_p(j:j+2, 2);
                    v3_p = faces_p(j:j+2, 3);
                    % 计算面的法向量
                    normal_p = cross(v2_p - v1_p, v3_p - v2_p);
                    normal_p = normal_p / norm(normal_p); % 单位法向量
                    % 计算点到面的向量
                    vector_p = candidateCoords - v1_p';
                    % 计算点到面的向量与法向量的点积
                    dotProduct_p = dot(vector_p, normal_p');
                    % 如果点在面的内部，则增加内部计数
                    if dotProduct_p <= 0
                        countInside_p = countInside_p + 1;
                    end
                end
                % 如果点在所有面的内部，标记为无效
                if countInside_p == 5
                    isValid = false;
                    break; % 退出当前四棱锥的检测
                end
            end
        end
    end
    % toc;

     if isValid
        % 检查长方体
        % tic;
        % for k = 1:numjinfei
            % 获取长方体的8个顶点
            vertices = jinfei(1, :);
            % 将 vertices 重新排列成 8x3 的矩阵
            vertices_matrix = reshape(vertices, 3, 8);
           
            
            % 计算长方体的距离
            distToBuilding = min(vecnorm(vertices_matrix - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToBuilding <= distanceThreshold
                 % 定义长方体的面
            faces = [
                vertices_matrix(:,[3, 2, 1, 4]); % 第一个面的四个顶点
                vertices_matrix(:,[5, 6, 7, 8]);  % 第二个面的四个顶点
                vertices_matrix(:,[5, 1, 2, 6]);  % 第三个面的四个顶点
                vertices_matrix(:,[6, 2, 3, 7]);  % 第四个面的四个顶点
                vertices_matrix(:,[7, 3, 4, 8]);  % 第五个面的四个顶点
                vertices_matrix(:,[8, 4, 1, 5]);  % 第六个面的四个顶点
            ];
                countInside = 0;
                for j = 1:3:18
                    % 面的三个顶点
                    v1 = faces(j:j+2, 1);
                    v2 = faces(j:j+2, 2);
                    v3 = faces(j:j+2, 3);
                    % 计算面的法向量
                    normal = cross(v2 - v1, v3 - v2);
                    normal = normal / norm(normal); % 单位法向量
                    % 计算点到面的向量
                    vector = candidateCoords - v1';
                    % 计算点到面的向量与法向量的点积
                    dotProduct = dot(vector, normal');
                    % 如果点在面的内部，则增加内部计数
                    if dotProduct <= 0
                        countInside = countInside + 1;
                    end
                end
                % 如果点在所有六个面的内部，标记为无效
                if countInside == 6
                    isValid = false;
                    % break; % 退出当前长方体的检测
                end
            end
        % end
     end
     % toc;
    if isValid
        % 检查长方体
        % tic;
        for k = 1:numBuildings
            % 获取长方体的8个顶点
            vertices = buildings(k, :);
            % 将 vertices 重新排列成 8x3 的矩阵
            vertices_matrix = reshape(vertices, 3, 8);
           
            
            % 计算长方体的距离
            distToBuilding = min(vecnorm(vertices_matrix - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToBuilding <= distanceThreshold
                 % 定义长方体的面
            faces = [
                vertices_matrix(:,[3, 2, 1, 4]); % 第一个面的四个顶点
                vertices_matrix(:,[5, 6, 7, 8]);  % 第二个面的四个顶点
                vertices_matrix(:,[5, 1, 2, 6]);  % 第三个面的四个顶点
                vertices_matrix(:,[6, 2, 3, 7]);  % 第四个面的四个顶点
                vertices_matrix(:,[7, 3, 4, 8]);  % 第五个面的四个顶点
                vertices_matrix(:,[8, 4, 1, 5]);  % 第六个面的四个顶点
            ];
                countInside = 0;
                for j = 1:3:18
                    % 面的三个顶点
                    v1 = faces(j:j+2, 1);
                    v2 = faces(j:j+2, 2);
                    v3 = faces(j:j+2, 3);
                    % 计算面的法向量
                    normal = cross(v2 - v1, v3 - v2);
                    normal = normal / norm(normal); % 单位法向量
                    % 计算点到面的向量
                    vector = candidateCoords - v1';
                    % 计算点到面的向量与法向量的点积
                    dotProduct = dot(vector, normal');
                    % 如果点在面的内部，则增加内部计数
                    if dotProduct <= 0
                        countInside = countInside + 1;
                    end
                end
                % 如果点在所有六个面的内部，标记为无效
                if countInside == 6
                    isValid = false;
                    break; % 退出当前长方体的检测
                end
            end
        end
    end
        if isValid
        validPoints(i) = true;
        end
    end
end
% end
% end
