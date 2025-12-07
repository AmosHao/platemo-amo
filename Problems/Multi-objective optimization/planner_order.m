classdef planner_order < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            obj.n=20;
            obj.m=7;
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = obj.n + 2 * (obj.m + 1); end
            obj.lower = ones(1,obj.D);
            obj.upper = obj.n*ones(1,obj.D);
            obj.encoding = 5*ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        M=obj.M;
        PopObj =zeros(N,M);
        %参数设置
        %直线距离
        % fl=[0.00, 6.96, 10.54, 10.31, 8.55, 4.23;
        %     6.96, 0.00, 3.90, 6.27, 7.90, 4.33;
        %     10.54, 3.90, 0.00, 4.70, 8.25, 6.98;
        %     10.31, 6.27, 4.70, 0.00, 4.26, 6.09;
        %     8.55, 7.90, 8.25, 4.26, 0.00, 5.00;
        %     4.23, 4.33, 6.98, 6.09, 5.00, 0.00
        %     ];
        %obj最小的距离
        fl=load('obj3min_f1.mat').f1;
        %无人机自身重量
        UAV_m=3;
        % 10个客户点的货物重量 
        demand = [3; 1.74; 0.2; 0.9; 2.0]; % 示例数据
        % demand = [3; 1.5; 1; 2; 1.5]; % 示例数据
        %
        a=25.5;
        b=0.2;c=0.6;d=0.2;
        P_ac=350;P_vc=250;P_dec=100;
        %约束
        maxlength=30;
        maxload=5;
        maxEC=500000;



if size(PopDec,1)>1
for j = 1:N 
Penalty1=0; Penalty2=0; Penalty3=0; 

%%%%%分段    
row = PopDec(j, :);
% 查找 -1 和 0 的位置
minus_one_positions = find(row == -1);
zero_positions = find(row == 0);
% 合并 -1 和 0 的位置并排序
all_positions = sort([minus_one_positions, zero_positions]);
% 初始化当前子段
current_segment = [];
segments = {};
% 遍历所有位置
for k = 1:length(all_positions)
    % 获取当前位置
    pos = all_positions(k); 
    % 如果当前位置是 -1，则将其加入当前段
    if row(pos) == -1
        current_segment = [current_segment, -1];
    end
    % 如果当前段不为空，则将其添加到结果中
    if ~isempty(current_segment)
        segments{end+1} = current_segment;
        current_segment = [];
    end
    % 如果当前位置是 0，则将其加入当前段
    if row(pos) == 0
        current_segment = [current_segment, 0];
    end
    
    % 如果当前位置不是最后一个位置，则将其后到下一个位置之间的元素加入当前段
    if k < length(all_positions)
        next_pos = all_positions(k+1);
        current_segment = [current_segment, row(pos+1:next_pos-1)];
    else
        % 如果当前位置是最后一个位置，则将其后到末尾的元素加入当前段
        current_segment = [current_segment, row(pos+1:end)];
    end
end
% 如果当前段不为空，则将其添加到结果中
if ~isempty(current_segment)
    segments{end+1} = current_segment;
end
     
%%计算三架飞机不同的距离
% 初始化结果向量
distances = zeros(1, 3);
% 遍历每个路径段
for i = 1:3
    segment = segments{i};
    total_distance = 0;
    
    % 计算路径段的总距离
    for k = 1:length(segment)-1
        if segment(k) ~= -1 && segment(k+1) ~= -1
            total_distance = total_distance + fl(segment(k)+1, segment(k+1)+1);
        end
    end
    
    % 如果路径段以 -1 结束，返回起点
    if segment(end)== -1
        total_distance = total_distance + fl(segment(end-1)+1, segment(1)+1);
    end
    if total_distance>maxlength
        Penalty1=Penalty1+total_distance;%约束每架无人机最大航程
    end
    % 存储总距离
    distances(i) = total_distance;
end
%目标一
PopObj(j,1)=sum(distances);       
 
EC=zeros(1,3);
for i = 1:3
    segment = segments{i};
    cargo_load = zeros( 1,length(segment) - 1);
    % 计算每架飞机每段之间的距离和载重
    cut_distance=[];
    for k = 1:length(segment) - 1
        % from = segment(k);
        to = segment(k + 1);
        if to == -1
            % 如果从起点出发或返回起点，载货量为 0
            cargo_load(k) = 0+UAV_m;
        else
            % 否则，计算从 from 到 to 之间的客户点货物需求量的和
            cargo_load(k) = UAV_m+sum(demand(segment(k+1:end-1)));
        end
    
         if segment(k) ~= -1 && segment(k+1) ~= -1
            cut_distance = [cut_distance ,fl(segment(k)+1, segment(k+1)+1)];
         end
    
    end
         % 如果路径段以 -1 结束，返回起点
         if segment(end)== -1
            cut_distance = [cut_distance , fl(segment(end-1)+1, segment(1)+1)];
         end

         if cargo_load(1)>maxload+UAV_m
             Penalty2=Penalty2+cargo_load(1);%约束每架无人机最大载重
         end
    % VC=zeros(1 ,length(segment) - 1);
    cut_distance=cut_distance.*1000;
   VC = cargo_load.^(3/4) .* a^(1/4) .* (cut_distance .* b .* P_ac).^(1/2) ...
    + P_vc .* c .* cut_distance .* (cargo_load .* a).^(-1/2) ...
    + cargo_load.^(3/4) .* a^(1/4) .* (cut_distance .* d .* P_dec).^(1/2);
 
   EC(i)=sum(VC);
   if EC(i)>maxEC
       Penalty3=Penalty3+EC(i)/10000;
   end
end
%目标2
PopObj(j,2)=sum(EC);

PopObj(j,1)=PopObj(j,1)+10^7*(Penalty1+Penalty2+Penalty3);
PopObj(j,2)=PopObj(j,2)+10^7*(Penalty1+Penalty2+Penalty3);
end

else
            j=1;
Penalty1=0; Penalty2=0; Penalty3=0; 

%%%%%分段    
row = PopDec(j, :);
% 查找 -1 和 0 的位置
minus_one_positions = find(row == -1);
zero_positions = find(row == 0);
% 合并 -1 和 0 的位置并排序
all_positions = sort([minus_one_positions, zero_positions]);
% 初始化当前子段
current_segment = [];
segments = {};
% 遍历所有位置
for k = 1:length(all_positions)
    % 获取当前位置
    pos = all_positions(k); 
    % 如果当前位置是 -1，则将其加入当前段
    if row(pos) == -1
        current_segment = [current_segment, -1];
    end
    % 如果当前段不为空，则将其添加到结果中
    if ~isempty(current_segment)
        segments{end+1} = current_segment;
        current_segment = [];
    end
    % 如果当前位置是 0，则将其加入当前段
    if row(pos) == 0
        current_segment = [current_segment, 0];
    end
    
    % 如果当前位置不是最后一个位置，则将其后到下一个位置之间的元素加入当前段
    if k < length(all_positions)
        next_pos = all_positions(k+1);
        current_segment = [current_segment, row(pos+1:next_pos-1)];
    else
        % 如果当前位置是最后一个位置，则将其后到末尾的元素加入当前段
        current_segment = [current_segment, row(pos+1:end)];
    end
end
% 如果当前段不为空，则将其添加到结果中
if ~isempty(current_segment)
    segments{end+1} = current_segment;
end
     
%%计算三架飞机不同的距离
% 初始化结果向量
distances = zeros(1, 3);
% 遍历每个路径段
for i = 1:3
    segment = segments{i};
    total_distance = 0;
    
    % 计算路径段的总距离
    for k = 1:length(segment)-1
        if segment(k) ~= -1 && segment(k+1) ~= -1
            total_distance = total_distance + fl(segment(k)+1, segment(k+1)+1);
        end
    end
    
    % 如果路径段以 -1 结束，返回起点
    if segment(end)== -1
        total_distance = total_distance + fl(segment(end-1)+1, segment(1)+1);
    end
    if total_distance>maxlength
        Penalty1=Penalty1+total_distance;%约束每架无人机最大航程
    end
    % 存储总距离
    distances(i) = total_distance;
end
%目标一
PopObj(j,1)=sum(distances);       
 
EC=zeros(1,3);
for i = 1:3
    segment = segments{i};
    cargo_load = zeros( 1,length(segment) - 1);
    % 计算每架飞机每段之间的距离和载重
    cut_distance=[];
    for k = 1:length(segment) - 1
        % from = segment(k);
        to = segment(k + 1);
        if to == -1
            % 如果从起点出发或返回起点，载货量为 0
            cargo_load(k) = 0+UAV_m;
        else
            % 否则，计算从 from 到 to 之间的客户点货物需求量的和
            cargo_load(k) = UAV_m+sum(demand(segment(k+1:end-1)));
        end
    
         if segment(k) ~= -1 && segment(k+1) ~= -1
            cut_distance = [cut_distance ,fl(segment(k)+1, segment(k+1)+1)];
         end
    
    end
         % 如果路径段以 -1 结束，返回起点
         if segment(end)== -1
            cut_distance = [cut_distance , fl(segment(end-1)+1, segment(1)+1)];
         end

         if cargo_load(1)>maxload+UAV_m
             Penalty2=Penalty2+cargo_load(1);%约束每架无人机最大载重
         end
    % VC=zeros(1 ,length(segment) - 1);
    cut_distance=cut_distance.*1000;
   VC = cargo_load.^(3/4) .* a^(1/4) .* (cut_distance .* b .* P_ac).^(1/2) ...
    + P_vc .* c .* cut_distance .* (cargo_load .* a).^(-1/2) ...
    + cargo_load.^(3/4) .* a^(1/4) .* (cut_distance .* d .* P_dec).^(1/2);
 
   EC(i)=sum(VC);
   if EC(i)>maxEC
       Penalty3=Penalty3+EC(i)/10000;
   end
end
%目标2
PopObj(j,2)=sum(EC);

PopObj(j,1)=PopObj(j,1)+10^7*(Penalty1+Penalty2+Penalty3);
PopObj(j,2)=PopObj(j,2)+10^7*(Penalty1+Penalty2+Penalty3);

end
        end


        %% Generate points on the Pareto front
        % function R = GetOptimum(obj,N)
        % % 
        %          R=load('SMEAGLT5.pf');
        %          % R=[20,1,0.0001];
        % % 
        % end

     
end
end
 function dist = distancePointToPolygon(px, py, polyX, polyY)
    % polyX, polyY 是平行四边形的顶点坐标
    % 使用 MATLAB 的内置方法来进行多边形测试
    % 计算到每个边的距离
    distances = [];
    for j = 1:length(polyX)
        % 得到当前边的两个顶点
        pt1 = [polyX(j), polyY(j)];
        pt2 = [polyX(mod(j, length(polyX))+1), polyY(mod(j, length(polyY))+1)];
        
        % 计算边的分量
        edge = pt2 - pt1;
        lengthEdge = norm(edge);
        
        if lengthEdge == 0
            % 如果边的长度为 0，使用点到顶点的距离
            distances(end + 1) = norm([px, py] - pt1);
        else
            % 计算当前点到边的投影
            projection = (dot([px, py] - pt1, edge) / lengthEdge^2) * edge + pt1;       
            % 计算距离
            distances(end + 1) = norm([px, py] - projection);
        end
    end
    
    % 返回最小距离
    dist = min(distances);
end
