function check2(obj,dec,validtrait)
% 文件名
filename = 'city_environment_data_new14.xlsx';

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


% 读取禁飞区数据
jinfei = xlsread(filename, 'jinfei');
numjinfei = size(jinfei, 1);

% 创建图形窗口
figure;
hold on;
view(3);
axis equal;
% 绘制地面长方形区域
xLimits = [6500 8000];
yLimits = [7000 8700];
zValue = 0;
fill3([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
      [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
      [zValue zValue zValue zValue], ...
      [1 1 0.8], 'FaceAlpha', 0.5); % 淡黄色
% 绘制地面长方形区域
xLimits = [9100 10000];
yLimits = [6300 8600];
zValue = 0;
fill3([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
      [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
      [zValue zValue zValue zValue], ...
      [1 1 0.8], 'FaceAlpha', 0.5); % 淡黄色
% 绘制地面长方形区域
xLimits = [2000 3100];
yLimits = [1200 2200];
zValue = 0;
fill3([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
      [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
      [zValue zValue zValue zValue], ...
      [1 1 0.8], 'FaceAlpha', 0.5); % 淡黄色
% 绘制地面长方形区域
xLimits = [6800 8100];
yLimits = [3300 4400];
zValue = 0;
fill3([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
      [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
      [zValue zValue zValue zValue], ...
      [1 1 0.8], 'FaceAlpha', 0.5); % 淡黄色
% 绘制地面长方形区域(river)
% 定义顶点坐标
x = [0, 0, 2000, 2000]; % x 坐标
y = [7300, 7450, 6050,5900]; % y 坐标
z = [0, 0, 0, 0]; % y 坐标
% 绘制平行四边形
fill3(x, y,z, [0.678, 0.847, 0.902], 'EdgeColor', 'k'); % 淡蓝色，边缘颜色为黑色

% 绘制地面圆形区域
circleCenter = [4600, 4700];
radius = 900;
theta = linspace(0, 2*pi, 100); % 生成圆形的角度数据
xCircle = circleCenter(1) + radius * cos(theta); % 生成圆形的 x 坐标
yCircle = circleCenter(2) + radius * sin(theta); % 生成圆形的 y 坐标
zCircle = zeros(size(theta)); % 圆形的 z 坐标为 0

fill3(xCircle, yCircle, zCircle, [0.6 1 0.6], 'FaceAlpha', 0.5); % 淡绿色
% 绘制长方体
for i = 1:numBuildings
    vertices = reshape(buildings(i, :), 3, 8)';
    vertices(:, 3) = vertices(:, 3) ;  % 将 z 坐标乘以 20
    % 确保填充的顺序正确
    fill3(vertices([1,2,6,5],1), vertices([1,2,6,5],2), vertices([1,2,6,5],3), [0.7 0.7 0.7], 'FaceAlpha', 0.8); % 灰色
    fill3(vertices([2,3,7,6],1), vertices([2,3,7,6],2), vertices([2,3,7,6],3), [0.7 0.7 0.7], 'FaceAlpha', 0.8); % 灰色
    fill3(vertices([3,4,8,7],1), vertices([3,4,8,7],2), vertices([3,4,8,7],3), [0.7 0.7 0.7], 'FaceAlpha', 0.8); % 灰色
    fill3(vertices([4,1,5,8],1), vertices([4,1,5,8],2), vertices([4,1,5,8],3), [0.7 0.7 0.7], 'FaceAlpha', 0.8); % 灰色
    fill3(vertices([1,2,3,4],1), vertices([1,2,3,4],2), vertices([1,2,3,4],3), [0.7 0.7 0.7], 'FaceAlpha', 0.8); % 灰色
end
% 绘制禁飞区
for i = 1:numjinfei
    vertices = reshape(jinfei(i, :), 3, 8)';
    vertices(:, 3) = vertices(:, 3);  % 保留原始 z 坐标
    % 确保填充的顺序正确
    fill3(vertices([1,2,6,5],1), vertices([1,2,6,5],2), vertices([1,2,6,5],3), [1 0 0], 'FaceAlpha', 0.8);  % 红色
    fill3(vertices([2,3,7,6],1), vertices([2,3,7,6],2), vertices([2,3,7,6],3), [1 0 0],'FaceAlpha', 0.8);  % 红色
    fill3(vertices([3,4,8,7],1), vertices([3,4,8,7],2), vertices([3,4,8,7],3), [1 0 0], 'FaceAlpha', 0.8);  % 红色
    fill3(vertices([4,1,5,8],1), vertices([4,1,5,8],2), vertices([4,1,5,8],3), [1 0 0], 'FaceAlpha', 0.8);  % 红色
    fill3(vertices([1,2,3,4],1), vertices([1,2,3,4],2), vertices([1,2,3,4],3), [1 0 0], 'FaceAlpha', 0.8);  % 红色
end
% 绘制圆柱体
for i = 1:numCylinders
    % 获取圆柱体的基本数据
    radius = cylinders(i, 3);  % 圆柱体的半径
    height = cylinders(i, 4);  % 圆柱体的高度
    x_center = cylinders(i, 1); % 圆心的 x 坐标
    y_center = cylinders(i, 2); % 圆心的 y 坐标
    
    % 生成圆柱体的侧面数据
    [X, Y, Z] = cylinder(radius, 20);
    % 调整 X 和 Y 坐标，使得 Y 为 0 的列中，X 对应为半径，Y 为半径的列中，X 对应为 0
    X(Y == 0) = radius; % Y 为 0 的列，X 对应为半径
    Y(X == 0) = radius; % X 为 0 的列，Y 对应为半径
    Z = Z * height; % 调整 Z 高度
    X = X + x_center; % 平移 X 坐标
    Y = Y + y_center; % 平移 Y 坐标
    
    % 绘制圆柱体的侧面
    surf(X, Y, Z, 'FaceColor', 'b', 'FaceAlpha', 0.8); % 蓝色侧面
    % % 生成顶面的数据
    % top_Z = height * ones(size(X)); % 顶面的 Z 坐标
    % 
    % % 绘制顶面
    % surf(X, Y, top_Z, 'FaceColor', 'r', 'FaceAlpha', 0.8); % 红色顶面

    % % 生成圆顶面的数据
    % theta = linspace(0, 2*pi, 100); % 角度范围
    % x_circle = radius * cos(theta) + x_center; % 圆的 X 坐标
    % y_circle = radius * sin(theta) + y_center; % 圆的 Y 坐标
    % z_circle = height * ones(size(x_circle)); % 顶面的 Z 坐标

    % 绘制圆柱体的顶面
    patch(X(2,:), Y(2,:), Z(2,:), 'b', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 蓝色顶面，去除边框

    % % （可选）绘制圆柱体的底面
    % bottom_x_circle = x_circle;
    % bottom_y_circle = y_circle;
    % bottom_z_circle = zeros(size(x_circle)); % 底面的 Z 坐标
    % patch(bottom_x_circle, bottom_y_circle, bottom_z_circle, 'b', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 蓝色底面，去除边框
end


% 绘制四棱锥
for i = 1:numPyramids
    vertices = reshape(pyramids(i, :), 3, 5)';
    vertices(:, 3) = vertices(:, 3) ;  % 将 z 坐标乘以 20
    fill3(vertices([1,2,3,4],1), vertices([1,2,3,4],2), vertices([1,2,3,4],3), [1 0.5 0], 'FaceAlpha', 0.8); % 橙色
    fill3([vertices(1,1), vertices(2,1), vertices(5,1)], [vertices(1,2), vertices(2,2), vertices(5,2)], [vertices(1,3), vertices(2,3), vertices(5,3)], [1 0.5 0], 'FaceAlpha', 0.8); % 橙色
    fill3([vertices(2,1), vertices(3,1), vertices(5,1)], [vertices(2,2), vertices(3,2), vertices(5,2)], [vertices(2,3), vertices(3,3), vertices(5,3)], [1 0.5 0], 'FaceAlpha', 0.8); % 橙色
    fill3([vertices(3,1), vertices(4,1), vertices(5,1)], [vertices(3,2), vertices(4,2), vertices(5,2)], [vertices(3,3), vertices(4,3), vertices(5,3)], [1 0.5 0], 'FaceAlpha', 0.8); % 橙色
    fill3([vertices(4,1), vertices(1,1), vertices(5,1)], [vertices(4,2), vertices(1,2), vertices(5,2)], [vertices(4,3), vertices(1,3), vertices(5,3)], [1 0.5 0], 'FaceAlpha', 0.8); % 橙色
end

% 绘制球体
for i = 1:numSpheres
    [X, Y, Z] = sphere(20);  % 生成单位球体网格数据
    X = X * spheres(i, 4) + spheres(i, 1);
    Y = Y * spheres(i, 4) + spheres(i, 2);
    Z = Z * spheres(i, 4) + spheres(i, 3);  % 球体底部在 z=0 以上
    
    % 绘制球体
    surf(X, Y, Z, 'FaceColor', 'r', 'FaceAlpha', 0.8); % 红色
end
% 绘制紫色小圆点
points = [500 500; 7000 3000; 9500 6000; 6050 9200; 1800 8950; 2800 4050];
scatter3(points(:,1), points(:,2), zeros(size(points, 1), 1), 60, [0.8, 0.6, 0.8], 'filled', 'MarkerEdgeColor', 'k'); % 紫色小圆点

% % 文件名
% filename = 'travel06.xlsx';
% obj = xlsread(filename, 'obj_data');
% dec = xlsread(filename, 'dec_data');
% [minValue, minIndex] = min(obj(:, 1));%obj第一列最小值
% dec1=[0.2,dec(minIndex,6:23),0.2,0.2,dec(minIndex,24:41),0.2,0.8,dec(minIndex,42:59),0.8];%dec第1行
% h1=plot3(dec1(:,1:20),dec1(:,21:40),dec1(:,41:60),'bo','MarkerSize',4,'MarkerFaceColor','g');
% hold on;
% %画线
%         % 坐标映射表 (坐标点与编号的关系)
% coordinates = [2.5, 0.5, 0.8, 0.2; 
%                3.5, 3, 0.8, 0.5; 
%                2.5, 5, 0.8, 0.8; 
%                2.5, 6.8, 0.8, 0.3; 
%                0.5, 5, 0.8, 1.0];  % 5个客户点的坐标和货物重量
%         indices = dec(minIndex, 1:5); % 前五列的编号
%         % 提取每个编号对应的坐标
%         coords = zeros(5, 3); % 5个点，每个点有3个坐标
%         for k = 1:5
%             coords(k, :) = coordinates(indices(k), 1:3);
%         end
%         xStart=[0.2,0.2,0.8];
%         z=dec(minIndex,42:59); y=dec(minIndex,24:41); x=dec(minIndex,6:23);
%         x=[xStart(1),x(1:3),coords(1,1),x(4:6),coords(2,1),x(7:9),coords(3,1),x(10:12),coords(4,1),x(13:15),coords(5,1),x(16:18),xStart(1)];
%         y=[xStart(2),y(1:3),coords(1,2),y(4:6),coords(2,2),y(7:9),coords(3,2),y(10:12),coords(4,2),y(13:15),coords(5,2),y(16:18),xStart(2)];
%         z=[xStart(3),z(1:3),coords(1,3),z(4:6),coords(2,3),z(7:9),coords(3,3),z(10:12),coords(4,3),z(13:15),coords(5,3),z(16:18),xStart(3)];
% 
% plot3(x(1:5), y(1:5), z(1:5), '-o', 'LineWidth',1, 'MarkerSize', 4, 'MarkerFaceColor', 'k','Color', 'r');
% plot3(x(5:9), y(5:9), z(5:9), '-o', 'LineWidth',1, 'MarkerSize', 4, 'MarkerFaceColor', 'k','Color', 'y');
% plot3(x(9:13), y(9:13), z(9:13), '-o', 'LineWidth',1, 'MarkerSize', 4, 'MarkerFaceColor', 'k','Color', 'g');
% plot3(x(13:17), y(13:17), z(13:17), '-o', 'LineWidth',1, 'MarkerSize', 4, 'MarkerFaceColor', 'k','Color', 'b');
% plot3(x(17:21), y(17:21), z(17:21), '-o', 'LineWidth',1, 'MarkerSize', 4, 'MarkerFaceColor', 'k','Color', 'm');
% plot3(x(21:25), y(21:25), z(21:25), '-o', 'LineWidth',1, 'MarkerSize', 4, 'MarkerFaceColor', 'k','Color', 'c');

% 文件名
% n_dots=165/3;
% filename = 'point2_penal2_dmin0.3_travel_singlepath1015_F3_penaly0_guideclose_maxfe20000_initialy_rand0.3.xlsx';
 % filename = 'maxhjddot2_disselect_initialselect_objduiying01.xlsx';

% filename = 'bj_dot1_1_dot2_6.xlsx';
% obj = xlsread(filename, 'obj_data');
% dec = xlsread(filename, 'dec_data');
% validtrait = xlsread(filename, 'validtrait_data');
% bj = xlsread(filename, 'bj_data');
vardim=size(dec,2);
n_dots=vardim/3;

% pop=dec;
% figure;
% hold on;
% vardim=size(dec,2);
% n_dots=vardim/3;
% view(3);
% axis equal;
% % for i=1:size(pop,1)
% dec1=pop(:,1:n_dots);
% dec2=pop(:,n_dots+1:2*n_dots);
% dec3=pop(:,2*n_dots+1:3*n_dots);
% plot3(dec1,dec2,dec3,'bo','MarkerSize',3,'MarkerFaceColor','g');
% hold on;
% % end
% hold on;
obj_valid=obj(find(validtrait(:,1)==1),:);
dec_valid=dec(find(validtrait(:,1)==1),:);
% bj_valid=bj(find(validtrait(:,1)==1),:);
[minValue, minIndex] = min(obj_valid(:, 3));%obj第一列最小值
% [minValue, minIndex] = min(obj(:, 1));%obj第一列最小值
% minIndex=6;
% dec1=[500,dec_valid(minIndex,1:n_dots),9500];
% dec2=[500,dec_valid(minIndex,n_dots+1:2*n_dots),6000];
% dec3=[800,dec_valid(minIndex,2*n_dots+1:3*n_dots),800];
% dec1=[500,dec(minIndex,1:n_dots),1800];
% dec2=[500,dec(minIndex,n_dots+1:2*n_dots),8900];
% dec3=[800,dec(minIndex,2*n_dots+1:3*n_dots),800];
dec1=[dec_valid(minIndex,1:n_dots)];
dec2=[dec_valid(minIndex,n_dots+1:2*n_dots)];
dec3=[dec_valid(minIndex,2*n_dots+1:3*n_dots)];
% dec1=dec1(:,find(bj_valid(minIndex,:)==0));
% dec2=dec2(:,find(bj_valid(minIndex,:)==0));
% dec3=dec3(:,find(bj_valid(minIndex,:)==0));
% dec1=[dec(minIndex,1:n_dots)];
% dec2=[dec(minIndex,n_dots+1:2*n_dots)];
% dec3=[dec(minIndex,2*n_dots+1:3*n_dots)];
%h1=plot3(dec1,dec2,dec3,'bo','MarkerSize',4,'MarkerFaceColor','g');
plot3(dec1,dec2,dec3, '-o', 'LineWidth',0.5, 'MarkerSize', 2, 'MarkerFaceColor', 'r','Color', 'g');
hold on;
% 画线
% 设置坐标轴属性

xlabel('X');
ylabel('Y');
zlabel('Z');
title('城市环境模型');

hold off;
end