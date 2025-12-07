function [validPoints]=checknode_singlepath_closefast(currentSol,varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids)
% %读取cubes数据
%  % 文件名
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

n_dots=varDim/3;
x=currentSol(:,1:n_dots);
y=currentSol(:,n_dots+1:2*n_dots);
z=currentSol(:,2*n_dots+1:3*n_dots);

validPoints = false(n_dots, 1);
distanceThreshold = 1;  % 例如：5单位

for i = 1:n_dots
    candidateCoords = [x(:, i), y(:, i), z(:, i)];
    isValid = true; % 默认认为是有效（在外部）

    % 检查球体
    tic;  % 开始计时
    for k = 1:numSpheres
        center = spheres(k, 1:3);
        radius = spheres(k, 4);
        distToCenter = norm(candidateCoords - center);
        if distToCenter <= radius
            isValid = false;
            break; % 退出当前球体的检测
        end
    end
    toc;  % 结束计时，输出时间

    if isValid
        % 检查圆柱体
        tic;  % 开始计时
        for k = 1:numCylinders
            center = cylinders(k, 1:2);
            radius = cylinders(k, 3);
            height = cylinders(k, 4);
            distToCenter = norm(candidateCoords(1:2) - center(1:2));
            if distToCenter <= radius && candidateCoords(3) <= height
                isValid = false;
                break; % 退出当前圆柱体的检测
            end
        end
        toc;  % 结束计时

    end
    
    if isValid
        % 检查四棱锥
        tic;  % 开始计时
        for k = 1:numPyramids
            vertices_p = pyramids(k, :);
            vertices_p = reshape(vertices_p, 3, 5);
            distToPyramid = min(vecnorm(vertices_p - candidateCoords', 2, 1));
            if distToPyramid <= distanceThreshold
                faces_p = [
                    vertices_p(:,[3, 2, 1]); 
                    vertices_p(:,[1,5,4]);  
                    vertices_p(:,[4,5,3]);  
                    vertices_p(:,[2,5,1]);  
                    vertices_p(:,[3,5,2]);  
                ];
                countInside_p = 0;
                for j = 1:3:15
                    v1_p = faces_p(j:j+2, 1);
                    v2_p = faces_p(j:j+2, 2);
                    v3_p = faces_p(j:j+2, 3);
                    normal_p = cross(v2_p - v1_p, v3_p - v2_p);
                    normal_p = normal_p / norm(normal_p);
                    vector_p = candidateCoords - v1_p';
                    dotProduct_p = dot(vector_p, normal_p');
                    if dotProduct_p <= 0
                        countInside_p = countInside_p + 1;
                    end
                end
                if countInside_p == 5
                    isValid = false;
                    break; 
                end
            end
        end
        toc;  % 结束计时
    end

    if isValid
        % 检查长方体
        tic;  % 开始计时
        for k = 1:numBuildings
            vertices = buildings(k, :);
            vertices_matrix = reshape(vertices, 3, 8);
            distToBuilding = min(vecnorm(vertices_matrix - candidateCoords', 2, 1));
            if distToBuilding <= distanceThreshold
                faces = [
                    vertices_matrix(:,[3, 2, 1, 4]); 
                    vertices_matrix(:,[5, 6, 7, 8]);  
                    vertices_matrix(:,[5, 1, 2, 6]);  
                    vertices_matrix(:,[6, 2, 3, 7]);  
                    vertices_matrix(:,[7, 3, 4, 8]);  
                    vertices_matrix(:,[8, 4, 1, 5]);  
                ];
                countInside = 0;
                for j = 1:3:18
                    v1 = faces(j:j+2, 1);
                    v2 = faces(j:j+2, 2);
                    v3 = faces(j:j+2, 3);
                    normal = cross(v2 - v1, v3 - v2);
                    normal = normal / norm(normal);
                    vector = candidateCoords - v1';
                    dotProduct = dot(vector, normal');
                    if dotProduct <= 0
                        countInside = countInside + 1;
                    end
                end
                if countInside == 6
                    isValid = false;
                    break; 
                end
            end
        end
        toc;  % 结束计时
    end

    if isValid
        validPoints(i) = true;
    end
end