function DIS = checknode_singlepath_close_dis(currentSol,varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei)
    
n_dots=varDim/3;
x=currentSol(:,1:n_dots);
y=currentSol(:,n_dots+1:2*n_dots);
z=currentSol(:,2*n_dots+1:3*n_dots);
dsafe=250;
DIS=0;
% 设置距离阈值（可以根据需要进行调整）
distanceThreshold = 1400;  % 例如：5单位
for i = 1:n_dots
    % 当前候选点坐标
    candidateCoords = [x(:, i), y(:, i), z(:, i)];
    for k = 1:numSpheres
        % 获取球体参数
        center = spheres(k, 1:3);
        % 计算点到球体圆心的距离
        distToCenter = norm(candidateCoords - center);
        dis = max(0, dsafe - distToCenter); 
        DIS = DIS + dis;
    end
    for k = 1:numCylinders
            % 获取圆柱体参数
            center = cylinders(k, 1:2);
            % 计算点到圆柱体圆心的距离
            distToCenter = norm(candidateCoords(1:2) - center(1:2));
            dis = max(0, dsafe - distToCenter); 
            DIS = DIS + dis;
    end
     for k = 1:numPyramids
            % 获取四棱锥的5个顶点
            vertices_p = pyramids(k, :);
            vertices_p = reshape(vertices_p, 3, 5);
             % 计算每个顶点的距离
            distToPyramid = min(vecnorm(vertices_p - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToPyramid <= distanceThreshold
            for j=1:5
                distToCenter = norm(candidateCoords - vertices_p(:,j));
            dis = max(0, dsafe - distToCenter); 
            DIS = DIS + dis;
            end
            end
     end
     for k = 1:numBuildings
            % 获取长方体的8个顶点
            vertices = buildings(k, :);
            % 将 vertices 重新排列成 8x3 的矩阵
            vertices_matrix = reshape(vertices, 3, 8);
             % 计算长方体的距离
            distToBuilding = min(vecnorm(vertices_matrix - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToBuilding <= distanceThreshold
            for j=1:8
            distToCenter = norm(candidateCoords - vertices_matrix(:,j));
            dis = max(0, dsafe - distToCenter); 
            DIS = DIS + dis;
            end
            end
     end
     vertices = jinfei(1, :);
            % 将 vertices 重新排列成 8x3 的矩阵
            vertices_matrix = reshape(vertices, 3, 8);
             % 计算长方体的距离
            distToBuilding = min(vecnorm(vertices_matrix - candidateCoords', 2, 1));
            
            % 只有在距离较近的情况下检查碰撞
            if distToBuilding <= distanceThreshold
            for j=1:8
            distToCenter = norm(candidateCoords - vertices_matrix(:,j));
            dis = max(0, dsafe - distToCenter); 
            DIS = DIS + dis;
            end
            end
                     
end
