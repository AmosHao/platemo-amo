classdef planner_order_v3 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            obj.n=10;  % 客户点数量
            obj.m=3;   % 无人机数量
            if isempty(obj.M); obj.M = 2; end
            % 编码维度：2*n个点（n个商家+n个客户点）+ (m-1)个0分隔符
            if isempty(obj.D); obj.D = 2*obj.n + (obj.m - 1); end
            obj.lower = ones(1,obj.D);
            obj.upper = 2*obj.n*ones(1,obj.D);  % 最大编号为20（商家编号1-10，客户点编号11-20）
            obj.encoding = 5*ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;%决策变量维度
        M=obj.M;%目标维度
        n=obj.n;%客户点数量
        m=obj.m;%无人机数量
        PopObj =zeros(N,M);
        % 参数设置
        % 无人机自身重量
        UAV_m=2;
        % 10个客户点的货物重量 --最后补0代表返回配送中心时的载重
        % 注意：商家取餐时载重增加，客户点配送时载重减少
        demand_q = [0.5; 0.7; 0.65; 0.9; 0.35; 0.2; 0.9; 0.1; 0.45; 0.77; 0]; 
        % 水平飞行速度
        v_h=10;
        % 垂直上升速度
        v_u=5;
        % 垂直下降速度
        v_d=3;
        % 重力加速度
        g=9.8;
        % 垂直高度
        h_up_down=40;
        % 引入数据文件（包含坐标、障碍物、惩罚系数等）
        % 数据文件路径：从当前文件位置到forOrderNew26文件夹
        current_file_dir = fileparts(mfilename('fullpath'));
        % 从 PlatEMO/Problems/Multi-objective optimization/Amos/ 到 PlatEMO/ 目录
        platemo_root = fileparts(fileparts(fileparts(current_file_dir)));
        data_file_path = fullfile(platemo_root, 'forOrderNew26', 'order_data.m');
        
        if exist(data_file_path, 'file')
            run(data_file_path);
        else
            % 如果找不到数据文件，使用默认数据（向后兼容）
            warning('找不到数据文件 %s，使用默认数据', data_file_path);
            // dotss = [
            //     5000, 5000, 800;   % 配送中心（在正中间）
            //     8156, 1270, 800;   % 商家1
            //     9058, 2785, 800;   % 商家2
            //     8500, 1576, 800;   % 商家3
            //      975, 1419, 800;   % 商家4
            //     2785, 5469, 800;   % 商家5
            //     9575, 9649, 800;   % 商家6
            //     4854, 8003, 800;   % 商家7
            //     7922, 9595, 800;   % 商家8
            //      357, 6557, 800;   % 商家9
            //     8491, 9340, 800;   % 商家10
            //      500,  500, 800;   % 客户点11
            //     1500, 2800, 800;   % 客户点12
            //     7000, 3000, 800;   % 客户点13
            //     9000, 5500, 800;   % 客户点14
            //     1500, 4500, 800;   % 客户点15
            //     9500, 5500, 800;   % 客户点16
            //      500, 7000, 800;   % 客户点17
            //     6000, 6100, 800;   % 客户点18
            //     1800, 8950, 800;   % 客户点19
            //     6050, 9200, 800    % 客户点20
            // ];
            // forbidden_zones = [
            //     2000, 2000, 4000, 4000;
            //     6000, 1000, 8000, 3000;
            // ];
            // crowded_zones = [
            //     6500, 6500, 500;
            //     3000, 7000, 800;
            // ];
            // noise_zones = [
            //     4500, 3500, 5500, 4500;
            //     8500, 3500, 9500, 4500;
            // ];
            // penalty_forbidden = 2000;
            // penalty_crowded = 500;
            // penalty_noise = 1000;
        end
        
        % 计算考虑障碍物的距离矩阵（21×21）
        % 行/列索引：1=配送中心(0), 2-11=商家1-10, 12-21=客户点11-20
        d_matrix = calculate_distance_matrix_with_obstacles(...
            dotss, forbidden_zones, crowded_zones, noise_zones, ...
            penalty_forbidden, penalty_crowded, penalty_noise);
        % 时间窗矩阵（21行：配送中心，商家1-10，客户点11-20）
        % 行索引：1=配送中心, 2-11=商家1-10, 12-21=客户点11-20
        TW_sec = [
            0, 100000000000;   % 配送中心
             0*60, 15*60;      % 商家1（出餐时间窗）
             0*60, 10*60;     % 商家2
             0*60, 10*60;     % 商家3
             0*60, 10*60;     % 商家4
             0*60, 10*60;     % 商家5
             0*60, 15*60;     % 商家6
             0*60, 12*60;     % 商家7
             0*60, 10*60;     % 商家8
             0*60, 10*60;     % 商家9
             0*60, 10*60;     % 商家10
            20*60, 40*60;     % 客户点11
            10*60, 30*60;     % 客户点12
             5*60, 25*60;     % 客户点13
             0*60, 20*60;     % 客户点14
             5*60, 25*60;     % 客户点15
            20*60, 40*60;     % 客户点16
            15*60, 35*60;     % 客户点17
            10*60, 30*60;     % 客户点18
             0*60, 20*60;     % 客户点19
             5*60, 25*60      % 客户点20
        ];
        % 时间窗惩罚项
        lam1=0.5;  % 从1.5降低到0.5，减少提前到达的惩罚
        lam2=1.0;  % 从2.0降低到1.0，减少晚到达的惩罚
        % 航迹规划得到的航机长度矩阵
        % d_matrix=zeros(n+1,n+1);  % 注释掉这行，使用下面的欧氏距离矩阵
        % 使用欧氏距离矩阵作为基础距离
        % 如果需要航迹规划，可以在此基础上进行优化
        
        %约束
        maxload=3;%最大载荷
        maxEC=800000;%电池容量

    if size(PopDec,1)>1
        for j = 1:N 
        Penalty1=0; Penalty2=0; Penalty3=0;  % Penalty3用于顺序约束
        route=PopDec(j,:);
        
        % ---------- 0. 检查顺序约束：商家必须在对应客户点之前访问 ----------
        % 商家编号：1-10，对应客户点编号：11-20
        % 对于每个客户点i+10，如果路径中包含客户点i+10，则必须在其之前包含商家i
        for merchant_id = 1:n
            customer_id = merchant_id + n;  % 客户点编号
            merchant_pos = find(route == merchant_id);
            customer_pos = find(route == customer_id);
            
            % 如果路径中包含客户点但不包含对应商家，或者商家在客户点之后，则违反约束
            if ~isempty(customer_pos)
                if isempty(merchant_pos)
                    % 有客户点但没有对应商家
                    Penalty3 = Penalty3 + 1000;
                else
                    % 检查商家是否在客户点之前（考虑所有0分隔符分段）
                    % 找到客户点和商家所在的段
                    idx0 = find(route == 0);
                    idx0 = [0 idx0 numel(route)+1];
                    
                    customer_seg = 0;
                    merchant_seg = 0;
                    for seg_idx = 1:length(idx0)-1
                        seg_range = (idx0(seg_idx)+1):(idx0(seg_idx+1)-1);
                        if any(route(seg_range) == customer_id)
                            customer_seg = seg_idx;
                        end
                        if any(route(seg_range) == merchant_id)
                            merchant_seg = seg_idx;
                        end
                    end
                    
                    % 如果商家和客户点在同一段，检查顺序
                    if merchant_seg == customer_seg && merchant_seg > 0
                        seg_range = (idx0(merchant_seg)+1):(idx0(merchant_seg+1)-1);
                        merchant_pos_in_seg = find(route(seg_range) == merchant_id, 1);
                        customer_pos_in_seg = find(route(seg_range) == customer_id, 1);
                        if merchant_pos_in_seg >= customer_pos_in_seg
                            Penalty3 = Penalty3 + 1000;  % 商家在客户点之后，违反约束
                        end
                    elseif merchant_seg > customer_seg && customer_seg > 0
                        Penalty3 = Penalty3 + 1000;  % 商家段在客户点段之后，违反约束
                    end
                end
            end
        end
        
        % ---------- 1. 根据 0 切分m段 ----------
        idx0   = find(route == 0);
        idx0   = [0 idx0 numel(route)+1];          % 方便循环
        routes = cell(m,1);
        for k = 1:m
            routes{k} = route(idx0(k)+1 : idx0(k+1)-1);
        end
        
        % 路径按 0 分段后直接用于后续计算

        % --- 2. 初始化输出 ---
        T_all = zeros(1,m);  % 每架无人机总时间
        E_all = zeros(1,m);  % 每架无人机总能耗
        load_all = zeros(1,m); % 每架无人机初始载重
        
        % --- 3. 每架无人机逐点计算 ---
        for k = 1:m
            seg = routes{k}; % 路径段
            if isempty(seg)
                continue;
            end
        
            % 路径点：depot -> seg(1) -> seg(2) -> ... -> seg(end) -> depot
            npts = length(seg) + 2; % 段内点数+2（depot起点和终点）
        
            % 当前剩余载荷计算（考虑商家取餐和客户点配送）
            % 商家（1-10）：取餐，载重增加
            % 客户点（11-20）：配送，载重减少
            segforload = [seg, 2*n+1];  % 添加配送中心作为终点（索引为2*n+1=21）
            load_kg = zeros(npts-1, 1);
            current_load = 0;  % 从配送中心出发时载重为0
            
            for i = 1:npts-1
                point_id = segforload(i);
                if point_id >= 1 && point_id <= n
                    % 商家（编号1-10）：取餐，载重增加
                    customer_idx = point_id;  % 商家i对应客户点i+10，但demand_q索引是1-10
                    current_load = current_load + demand_q(customer_idx);
                elseif point_id > n && point_id <= 2*n
                    % 客户点（编号11-20）：配送，载重减少
                    customer_idx = point_id - n;  % 客户点编号转换为demand_q索引（1-10）
                    current_load = current_load - demand_q(customer_idx);
                end
                load_kg(i) = current_load;
            end
            load_all(k) = max(load_kg);  % 保存每架无人机最大载重，以计算约束
        
               T_k_all = 0; 
               E_k     = 0;
               segforlength = [0,seg,0];
               % 累积时间（到达每个点的时间）
               time_cum = 0;
            
            for i = 1:npts-1 
                % 此处引入航迹规划算法，输入dot索引和坐标, 输出路程。需要缓存判断。
                % [length]=runPlatemo_forOrder_PRMO(problem,algorithm,N,maxFE,s,M,dots,dot1,dot2,whichObj);
                dot1=segforlength(i)+1;
                dot2=segforlength(i+1)+1; 
                % dots=[dotss(dot1,:),dotss(dot2,:),];
                % if d_matrix(dot1,dot2)==0
                % [f_length]=runPlatemo_forOrder_PRMO(@planner_simple_maxhjd_newobj,@PRMO,100,300,20,3,dots,dot1,dot2,1);
                % d_matrix(dot1,dot2)=f_length;
                % d_matrix(dot2,dot1)=f_length;
                % d_hor=f_length;
                % else
                d_hor=d_matrix(dot1,dot2);
                % end
                

                % 上升（仅在从 depot 或客户点起飞时）
                t_up = h_up_down / v_u;
                % 此处写入上升功率计算公式
                UAV_w = (UAV_m+load_kg(i))*g; %无人机自重w（T）
                v_0 = UAV_w/1.23235;
                P_up_i = 79.85628 + UAV_w*(v_u/2+sqrt((v_u/2)^2+v_0));
                E_up = P_up_i * t_up;
        
                % 水平飞行
                t_hor = d_hor / v_h;%此处d_hor为航迹规划算法输出的值或矩阵内存调用
                %此处写入水平功率计算公式
                P_hor_i = 79.85628*(1+0.00020833*v_h^4) + (1.1*UAV_w^(3/2)/sqrt(1.23235))*(sqrt((1+v_h^4)/(4*v_0^4))-v_h^2/(2*v_0))^(1/2) + 0.009249*v_h;
                E_hor = P_hor_i * t_hor;
        
                % 下降（仅在到达客户点或 depot 时）
                t_down = h_up_down / v_d;
                %此处写入下降功率计算公式
                v_d_0=v_d/v_0;
                P_down_i = 79.85625 + UAV_w*v_0*(0.974 -1.125*v_d_0 -1.372*v_d_0^2 -1.718*v_d_0^3 -0.655*v_d_0^4);
                E_down = P_down_i * t_down;

               % 计算去往每个点路上所花费的时间（单段）
               T_k_single = t_up + t_hor + t_down;
               % 累积到达时间（从出发到当前“到达点”的总时间）
               time_cum = time_cum + T_k_single;

               % ---- 基于"到达点"的时间窗计算时间成本 ----
               % 当前段是 segforlength(i) -> segforlength(i+1)
               % 到达点索引：0=配送中心，1-10=商家，11-20=客户点
               dest_idx = segforlength(i+1);
               if dest_idx == 0
                   % 配送中心对应 TW_sec 的第1行
                   tw_row = 1;
               elseif dest_idx >= 1 && dest_idx <= n
                   % 商家（1-10）对应 TW_sec 的第2-11行
                   tw_row = dest_idx + 1;
               elseif dest_idx > n && dest_idx <= 2*n
                   % 客户点（11-20）对应 TW_sec 的第12-21行
                   tw_row = dest_idx + 1;  % 客户点11对应第12行，客户点20对应第21行
               else
                   tw_row = 1;  % 默认使用配送中心时间窗
               end
               aa = TW_sec(tw_row,1);   % 时间窗下界
               bb = TW_sec(tw_row,2);   % 时间窗上界

               if time_cum < aa         % 提前到达
                   T_k = lam1*(aa - time_cum);
               elseif time_cum <= bb    % 在时间窗内到达
                   T_k = time_cum - aa;
               else                     % 晚到达
                   T_k = lam2*(time_cum - bb);
               end
        
               % 累加
                T_k_all = T_k_all + T_k;
                E_k = E_k + E_up + E_hor + E_down;
            end
        
            T_all(k) = T_k_all;%每架飞机的
            E_all(k) = E_k;
        end
        
        % 目标一
        PopObj(j,1)=sum(E_all); 
        % 目标二
        PopObj(j,2)=sum(T_all); 
         
        % 约束--载重
        for kk = 1:m
            load_k = load_all(kk);
            if load_k > maxload
                Penalty1=Penalty1+(load_k-maxload);
            end
        end
        
        % 约束--电量
        for kk = 1:m
            power_k = E_all(kk);
            if  power_k > maxEC
                Penalty2=Penalty2+((power_k-maxEC)/1000);%转为KJ，对应载重约束kg
            end
        end
        
        PopObj(j,1)=PopObj(j,1)+10*(Penalty1+Penalty2)+Penalty3;  % 添加顺序约束惩罚
        PopObj(j,2)=PopObj(j,2)+10*(Penalty1+Penalty2)+Penalty3;  % 添加顺序约束惩罚
        end
        
    else
                    j=1;
        Penalty1=0; Penalty2=0; Penalty3=0;  % Penalty3用于顺序约束
        route=PopDec(j,:);
        
        % ---------- 0. 检查顺序约束：商家必须在对应客户点之前访问 ----------
        % 商家编号：1-10，对应客户点编号：11-20
        for merchant_id = 1:n
            customer_id = merchant_id + n;  % 客户点编号
            merchant_pos = find(route == merchant_id);
            customer_pos = find(route == customer_id);
            
            if ~isempty(customer_pos)
                if isempty(merchant_pos)
                    Penalty3 = Penalty3 + 1000;
                else
                    idx0 = find(route == 0);
                    idx0 = [0 idx0 numel(route)+1];
                    
                    customer_seg = 0;
                    merchant_seg = 0;
                    for seg_idx = 1:length(idx0)-1
                        seg_range = (idx0(seg_idx)+1):(idx0(seg_idx+1)-1);
                        if any(route(seg_range) == customer_id)
                            customer_seg = seg_idx;
                        end
                        if any(route(seg_range) == merchant_id)
                            merchant_seg = seg_idx;
                        end
                    end
                    
                    if merchant_seg == customer_seg && merchant_seg > 0
                        seg_range = (idx0(merchant_seg)+1):(idx0(merchant_seg+1)-1);
                        merchant_pos_in_seg = find(route(seg_range) == merchant_id, 1);
                        customer_pos_in_seg = find(route(seg_range) == customer_id, 1);
                        if merchant_pos_in_seg >= customer_pos_in_seg
                            Penalty3 = Penalty3 + 1000;
                        end
                    elseif merchant_seg > customer_seg && customer_seg > 0
                        Penalty3 = Penalty3 + 1000;
                    end
                end
            end
        end
        
        % ---------- 1. 根据 0 切分m段 ----------
        idx0   = find(route == 0);
        idx0   = [0 idx0 numel(route)+1];          % 方便循环
        routes = cell(m,1);
        for k = 1:m
            routes{k} = route(idx0(k)+1 : idx0(k+1)-1);
        end
        
        % 路径按 0 分段后直接用于后续计算
        % --- 2. 初始化输出 ---
        T_all = zeros(1,m);  % 每架无人机总时间
        E_all = zeros(1,m);  % 每架无人机总能耗
        load_all = zeros(1,m); % 每架无人机初始载重
        
        % --- 3. 每架无人机逐点计算 ---
        for k = 1:m
            seg = routes{k}; % 1 2 3
            if isempty(seg)
                continue;
            end
        
            % 路径点：depot -> seg(1) -> seg(2) -> ... -> seg(end) -> depot
            % pts = [depot; customerXY(seg,:); depot];
            npts = length(seg) + 2; %3+2=5
        
            % 当前剩余载荷计算（考虑商家取餐和客户点配送）
            % 商家（1-10）：取餐，载重增加
            % 客户点（11-20）：配送，载重减少
            segforload = [seg, 2*n+1];  % 添加配送中心作为终点（索引为2*n+1=21）
            load_kg = zeros(npts-1, 1);
            current_load = 0;  % 从配送中心出发时载重为0
            
            for i = 1:npts-1
                point_id = segforload(i);
                if point_id >= 1 && point_id <= n
                    % 商家（编号1-10）：取餐，载重增加
                    customer_idx = point_id;  % 商家i对应客户点i+10，但demand_q索引是1-10
                    current_load = current_load + demand_q(customer_idx);
                elseif point_id > n && point_id <= 2*n
                    % 客户点（编号11-20）：配送，载重减少
                    customer_idx = point_id - n;  % 客户点编号转换为demand_q索引（1-10）
                    current_load = current_load - demand_q(customer_idx);
                end
                load_kg(i) = current_load;
            end
            load_all(k) = max(load_kg);  % 保存每架无人机最大载重，以计算约束
        
               T_k_all = 0; 
               E_k     = 0;
               segforlength = [0,seg,0];
               % 累积时间（到达每个点的时间）
               time_cum = 0;
            
            for i = 1:npts-1 
                % 此处引入航迹规划算法，输入dot索引和坐标, 输出路程。需要缓存判断。
                % [length]=runPlatemo_forOrder_PRMO(problem,algorithm,N,maxFE,s,M,dots,dot1,dot2,whichObj);
                dot1=segforlength(i)+1;
                dot2=segforlength(i+1)+1; 
                % dots=[dotss(dot1,:),dotss(dot2,:),];
                % if d_matrix(dot1,dot2)==0
                % [f_length]=runPlatemo_forOrder_PRMO(@planner_simple_maxhjd_newobj,@simplePRGOCEA,100,30000,20,3,dots,dot1,dot2,1);
                % d_matrix(dot1,dot2)=f_length;
                % d_matrix(dot2,dot1)=f_length;
                % d_hor=f_length;
                % else
                d_hor=d_matrix(dot1,dot2);
                % end
        
                % 上升（仅在从 depot 或客户点起飞时）
                t_up = h_up_down / v_u;
                % 此处写入上升功率计算公式
                UAV_w = (UAV_m+load_kg(i))*g; %无人机自重w（T）
                v_0 = UAV_w/1.23235;
                P_up_i = 79.85628 + UAV_w*(v_u/2+sqrt((v_u/2)^2+v_0));
                E_up = P_up_i * t_up;
        
                % 水平飞行
                t_hor = d_hor / v_h;%此处d_hor为航迹规划算法输出的值或矩阵内存调用
                %此处写入水平功率计算公式
                P_hor_i = 79.85628*(1+0.00020833*v_h^4) + (1.1*UAV_w^(3/2)/sqrt(1.23235))*(sqrt((1+v_h^4)/(4*v_0^4))-v_h^2/(2*v_0))^(1/2) + 0.009249*v_h;
                E_hor = P_hor_i * t_hor;
        
                % 下降（仅在到达客户点或 depot 时）
                t_down = h_up_down / v_d;
                %此处写入下降功率计算公式
                v_d_0=v_d/v_0;
                P_down_i = 79.85625 + UAV_w*v_0*(0.974 -1.125*v_d_0 -1.372*v_d_0^2 -1.718*v_d_0^3 -0.655*v_d_0^4);
                E_down = P_down_i * t_down;

               % 计算去往每个点路上所花费的时间（单段）
               T_k_single = t_up + t_hor + t_down;
               % 累积到达时间（从出发到当前“到达点”的总时间）
               time_cum = time_cum + T_k_single;

               % ---- 基于"到达点"的时间窗计算时间成本 ----
               % 当前段是 segforlength(i) -> segforlength(i+1)
               % 到达点索引：0=配送中心，1-10=商家，11-20=客户点
               dest_idx = segforlength(i+1);
               if dest_idx == 0
                   % 配送中心对应 TW_sec 的第1行
                   tw_row = 1;
               elseif dest_idx >= 1 && dest_idx <= n
                   % 商家（1-10）对应 TW_sec 的第2-11行
                   tw_row = dest_idx + 1;
               elseif dest_idx > n && dest_idx <= 2*n
                   % 客户点（11-20）对应 TW_sec 的第12-21行
                   tw_row = dest_idx + 1;  % 客户点11对应第12行，客户点20对应第21行
               else
                   tw_row = 1;  % 默认使用配送中心时间窗
               end
               aa = TW_sec(tw_row,1);   % 时间窗下界
               bb = TW_sec(tw_row,2);   % 时间窗上界

               if time_cum < aa         % 提前到达
                   T_k = lam1*(aa - time_cum);
               elseif time_cum <= bb    % 在时间窗内到达
                   T_k = time_cum - aa;
               else                     % 晚到达
                   T_k = lam2*(time_cum - bb);
               end
        
               % 累加
                T_k_all = T_k_all + T_k;
                E_k = E_k + E_up + E_hor + E_down;
            end
        
            T_all(k) = T_k_all;%每架飞机的
            E_all(k) = E_k;
        end
        
        % 目标一
        PopObj(j,1)=sum(E_all); 
        % 目标二
        PopObj(j,2)=sum(T_all); 
         
        % 约束--载重
        for kk = 1:m
            load_k = load_all(kk);
            if load_k > maxload
                Penalty1=Penalty1+(load_k-maxload);
            end
        end
        
        % 约束--电量
        for kk = 1:m
            power_k = E_all(kk);
            if  power_k > maxEC
                Penalty2=Penalty2+((power_k-maxEC)/1000);%转为KJ，对应载重约束kg
            end
        end
        
        % 调试输出：打印前几个解的目标值和惩罚值
        if j <= 3  % 只打印前3个解的信息
            fprintf('解 %d: E_all=[%.2f,%.2f,%.2f], T_all=[%.2f,%.2f,%.2f]\n', ...
                j, E_all(1), E_all(2), E_all(3), T_all(1), T_all(2), T_all(3));
            fprintf('解 %d: 原始目标值: E=%.2f, T=%.2f\n', j, PopObj(j,1), PopObj(j,2));
            fprintf('解 %d: 惩罚值: Penalty1=%.2f, Penalty2=%.2f, Penalty3=%.2f\n', j, Penalty1, Penalty2, Penalty3);
        end
        
        PopObj(j,1)=PopObj(j,1)+10*(Penalty1+Penalty2)+Penalty3;  % 添加顺序约束惩罚
        PopObj(j,2)=PopObj(j,2)+10*(Penalty1+Penalty2)+Penalty3;  % 添加顺序约束惩罚
        
        % 调试输出：打印最终目标值
        if j <= 3
            fprintf('解 %d: 最终目标值: E=%.2f, T=%.2f\n\n', j, PopObj(j,1), PopObj(j,2));
        end
        end
        end


        %% Generate points on the Pareto front
        % function R = GetOptimum(obj,N)
        % % 
        %          R=load('SMEAGLT5.pf');
        % % 
        % end

     
end
end

%% 测试函数：用于诊断目标函数计算问题（针对 v2 版本）
% function test_objective_function_v2()
%     % 创建一个简单的测试解
%     test_solution = [1, 2, 3, 0, 4, 5, 6, 0, 7, 8, 9, 10];
% 
%     % 创建问题实例
%     problem = planner_order_v2();
%     problem.Setting();
% 
%     % 计算目标值
%     obj = problem.CalObj(test_solution);
% 
%     fprintf('测试解: [');
%     fprintf('%d ', test_solution);
%     fprintf(']\n');
%     fprintf('目标值(v2): f1=%.2f, f2=%.2f\n', obj(1), obj(2));
% 
%     % 测试几个不同的解
%     test_solutions = [
%         [1, 2, 3, 0, 4, 5, 6, 0, 7, 8, 9, 10];
%         [10, 9, 8, 0, 7, 6, 5, 0, 4, 3, 2, 1];
%         [1, 3, 5, 0, 2, 4, 6, 0, 7, 9, 8, 10]
%     ];
% 
%     fprintf('\n测试多个解(v2):\n');
%     for i = 1:size(test_solutions, 1)
%         obj = problem.CalObj(test_solutions(i, :));
%         fprintf('解%d: f1=%.2f, f2=%.2f\n', i, obj(1), obj(2));
%     end
% end


