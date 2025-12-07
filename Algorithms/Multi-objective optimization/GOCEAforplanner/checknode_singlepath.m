function [validPoints]=checknode_singlepath(currentSol,varDim)
%读取cubes数据
 % 文件名
filename = 'city_environment_data_6.xlsx';
% 读取长方体数据
buildings = xlsread(filename, 'Buildings');
numBuildings = size(buildings, 1);
% 读取圆柱体数据
cylinders = xlsread(filename, 'Cylinders');
numCylinders = size(cylinders, 1);
% 读取四棱锥数据
pyramids = xlsread(filename, 'Pyramids');
numPyramids = size(pyramids, 1);
% 读取球体数据
spheres = xlsread(filename, 'Spheres');
numSpheres = size(spheres, 1);

% pop=xlsread('travel02', 'dec_data');
n_dots=varDim/3;
x=currentSol(:,1:n_dots);
y=currentSol(:,n_dots+1:2*n_dots);
z=currentSol(:,2*n_dots+1:3*n_dots);

% centerPoints = zeros(numBuildings, 3); % 存储长方体中心点的坐标
% for i = 1:numBuildings
%     % 获取当前长方体的8个顶点
%     vertices = reshape(buildings(i, :), [3, 8])';
%     % 计算中心点
%     centerPoints(i, :) = mean(vertices, 1);
% end
% % 计算每个航迹点到每个长方体中心点的距离
% numTrajectories = length(x); % 18
% distances = zeros(numTrajectories, numBuildings);
% for i = 1:numTrajectories
%     % 当前航迹点的坐标
%     trajectoryPoint = [x(i), y(i), z(i)];
%     % 计算当前航迹点到所有长方体中心点的距离
%     for j = 1:numBuildings
%         % 当前长方体中心点的坐标
%         centerPoint = centerPoints(j, :);
%         % 计算欧几里得距离
%         distances(i, j) = sqrt(sum((trajectoryPoint - centerPoint).^2));
%     end
% end

% for ii=1:size(pop,1)
% 进行碰撞检测
% 初始化标记数组，假设所有点初始为无效
validPoints = false(n_dots, 1);
% while any(~validPoints)
% for i = find(~validPoints)'
for i = 1:n_dots
    % 当前候选点坐标
    candidateCoords = [x(:,i), y(:,i), z(:,i)];
    % 检查该点是否与任何一个长方体发生碰撞
    isValid = false;
        % 检查球体
        for k = 1:numSpheres
            % 获取球体参数
            center = spheres(k, 1:3);
            radius = spheres(k, 4);
            
            % 计算点到球体圆心的距离
            distToCenter = norm(candidateCoords - center);
            % 检查点是否在球体内
            if distToCenter <= radius
                isValid = false;
                break; % 退出当前球体的检测
            else
                isValid = true;
            end
        end

    if isValid
        % 检查圆柱体
        for k = 1:numCylinders
            % 获取圆柱体参数
            center = cylinders(k, 1:2);
            radius = cylinders(k, 3);
            height = cylinders(k, 4);
            
            % 计算点到圆柱体圆心的距离
            distToCenter = norm(candidateCoords(1:2) - center(1:2));
            % 检查点是否在圆柱体的圆盘范围内
            if distToCenter <= radius && candidateCoords(3)<= height
                isValid = false;
                break; % 退出当前圆柱体的检测
            else
                isValid = true;
            end
        end
    end
    
    if isValid
        % 检查四棱锥
        for k = 1:numPyramids
        % 获取四棱锥的5个顶点
        vertices_p = pyramids(k, :);
        vertices_p = reshape(vertices_p, 3, 5);
            % 计算四棱锥的面
            % faces = [
            %     vertices([1, 2, 3], :); % 底面
            %     vertices([1, 2, 4], :); % 侧面1
            %     vertices([2, 3, 4], :); % 侧面2
            %     vertices([3, 1, 4], :); % 侧面3
            % ];
         faces_p= [
                vertices_p(:,[3, 2, 1]); % 底面
                vertices_p(:,[1,5,4]);  % 第二个面的四个顶点
                vertices_p(:,[4,5,3]);  % 第三个面的四个顶点
                vertices_p(:,[2,5,1]);  % 第四个面的四个顶点
                vertices_p(:,[3,5,2]);  % 第四个面的四个顶点
            ];
            % 检查点是否在四棱锥内
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
            else
                isValid = true;
            end
        end
    end

    if isValid
        % 检查长方体
    for k = 1:numBuildings
        % 获取长方体的8个顶点
        vertices = buildings(k, :);
        % 将 vertices 重新排列成 8x3 的矩阵
        vertices_matrix = reshape(vertices, 3, 8);
        % 定义长方体的面
        faces= [
                vertices_matrix(:,[3, 2, 1, 4]); % 第一个面的四个顶点
                vertices_matrix(:,[5, 6, 7, 8]);  % 第二个面的四个顶点
                vertices_matrix(:,[5, 1, 2, 6]);  % 第三个面的四个顶点
                vertices_matrix(:,[6, 2, 3, 7]);  % 第四个面的四个顶点
                vertices_matrix(:,[7, 3, 4, 8]);  % 第五个面的四个顶点
                vertices_matrix(:,[8, 4, 1, 5]);  % 第六个面的四个顶点
            ];
        % 检查点是否在长方体内
        % % 对每一个面进行检测
            countInside=0;
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
        else
            isValid = true;
        end
    end
    end
 % 如果点在所有形状的外部，标记为有效
        if isValid
            validPoints(i) = true;
        else
            % % 如果点仍然无效，对其进行随机扰动
            % perturbation = 0.1 * randn(1, 3); % 产生小的随机扰动
            % candidateCoords = candidateCoords + perturbation; % 添加扰动
            % 
            % % 将扰动后的坐标存回原数组中
            % x(ii,i) = candidateCoords(1);
            % y(ii,i) = candidateCoords(2);
            % z(ii,i) = candidateCoords(3);
        end
        validPoints=validPoints';
end
% end
% end
