function [currentSol,validPoints,bj]=checknode_singlepath_close_newbj_FANGZHEN(currentSol,varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2)
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
n_dots=varDim/3;
% x=[min(dot1(1),dot2(1)),currentSol(:,1:n_dots),max(dot1(1),dot2(1))];
% y=[min(dot1(2),dot2(2)),currentSol(:,n_dots+1:2*n_dots),max(dot1(2),dot2(2))];
% z=[min(dot1(3),dot2(3)),currentSol(:,2*n_dots+1:3*n_dots),max(dot1(3),dot2(3))];
x=[dot1(1),currentSol(:,1:n_dots),dot2(1)];
y=[dot1(2),currentSol(:,n_dots+1:2*n_dots),dot2(2)];
z=[dot1(3),currentSol(:,2*n_dots+1:3*n_dots),dot2(3)];


validPoints = false(n_dots, 1);
bj=false(n_dots, 1);
% 设置距离阈值（可以根据需要进行调整）
distanceThreshold = 1400;  % 例如：5单位

% AABB 预计算：用于“点在长方体内”的等价快速判定（包含高度 z）
xcols = 1:3:22; ycols = 2:3:23; zcols = 3:3:24;
if numBuildings > 0
    bxmin = min(buildings(:, xcols), [], 2); bxmax = max(buildings(:, xcols), [], 2);
    bymin = min(buildings(:, ycols), [], 2); bymax = max(buildings(:, ycols), [], 2);
    bzmin = min(buildings(:, zcols), [], 2); bzmax = max(buildings(:, zcols), [], 2);
else
    bxmin = []; bxmax = []; bymin = []; bymax = []; bzmin = []; bzmax = [];
end
if numjinfei > 0
    jxmin = min(jinfei(:, xcols), [], 2); jxmax = max(jinfei(:, xcols), [], 2);
    jymin = min(jinfei(:, ycols), [], 2); jymax = max(jinfei(:, ycols), [], 2);
    jzmin = min(jinfei(:, zcols), [], 2); jzmax = max(jinfei(:, zcols), [], 2);
else
    jxmin = []; jxmax = []; jymin = []; jymax = []; jzmin = []; jzmax = [];
end

for i = 1:n_dots
    % if i==7
    % aaa=113
    % end
    % 当前候选点坐标（第 i+1 个航迹点）
    candidateCoords = [x(:, i+1), y(:, i+1), z(:, i+1)];
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
            radius = cylinders(k, 3)*1.01;
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
        % 禁飞盒（jinfei）：AABB 等价点内判定（含 z）
        if ~isempty(jxmin)
            dx = max([jxmin - candidateCoords(1), zeros(size(jxmin)), candidateCoords(1) - jxmax], [], 2);
            dy = max([jymin - candidateCoords(2), zeros(size(jymin)), candidateCoords(2) - jymax], [], 2);
            dz = max([jzmin - candidateCoords(3), zeros(size(jzmin)), candidateCoords(3) - jzmax], [], 2);
            dmin = sqrt(dx.^2 + dy.^2 + dz.^2);
            close = (dmin <= distanceThreshold);
            if any(close)
                inside = close & ...
                    (candidateCoords(1) >= jxmin) & (candidateCoords(1) <= jxmax) & ...
                    (candidateCoords(2) >= jymin) & (candidateCoords(2) <= jymax) & ...
                    (candidateCoords(3) >= jzmin) & (candidateCoords(3) <= jzmax);
                if any(inside)
                    isValid = false;
                end
            end
        end
     end
     % toc;
    if isValid
        % 建筑长方体：AABB 等价点内判定（含 z），向量化一次判断所有楼
        if ~isempty(bxmin)
            dx = max([bxmin - candidateCoords(1), zeros(size(bxmin)), candidateCoords(1) - bxmax], [], 2);
            dy = max([bymin - candidateCoords(2), zeros(size(bymin)), candidateCoords(2) - bymax], [], 2);
            dz = max([bzmin - candidateCoords(3), zeros(size(bzmin)), candidateCoords(3) - bzmax], [], 2);
            dmin = sqrt(dx.^2 + dy.^2 + dz.^2);
            close = (dmin <= distanceThreshold);
            if any(close)
                inside = close & ...
                    (candidateCoords(1) >= bxmin) & (candidateCoords(1) <= bxmax) & ...
                    (candidateCoords(2) >= bymin) & (candidateCoords(2) <= bymax) & ...
                    (candidateCoords(3) >= bzmin) & (candidateCoords(3) <= bzmax);
                if any(inside)
                    isValid = false;
                end
            end
        end
    end
    % toc;

    % 如果点在所有形状的外部，标记为有效
    if isValid
        validPoints(i) = true;
        bj(i)=true;
        currentSol(:,i)=candidateCoords(1,1);
        currentSol(:,i+n_dots)=candidateCoords(1,2);
        currentSol(:,i+2*n_dots)=candidateCoords(1,3);
    end
end
% end
% end
