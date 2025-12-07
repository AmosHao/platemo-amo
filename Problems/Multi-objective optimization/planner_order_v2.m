classdef planner_order_v2 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            obj.n=10;
            obj.m=3;
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = obj.n + (obj.m - 1); end
            obj.lower = ones(1,obj.D);
            obj.upper = obj.n*ones(1,obj.D);
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
        % 客户点坐标矩阵
        % dotss = [
        %     5500, 5000, 800;   % 配送中心（第1行）
        %      500,  500, 800;   % 客户点1
        %     3800, 2800, 800;   % 客户点2
        %     7000, 3000, 800;   % 客户点3
        %     9000, 4000, 800;   % 客户点4
        %     2800, 4500, 800;   % 客户点5
        %     9500, 6000, 800;   % 客户点6
        %      500, 7000, 800;   % 客户点7
        %     6000, 6100, 800;   % 客户点8
        %     1800, 8950, 800;   % 客户点9
        %     6050, 9200, 800    % 客户点10
        % ];
        % 欧氏距离矩阵
        d_matrix = [
               0        6364.0   2651.9   3162.3   5000.0   2702.9   5025.0   6324.6   1581.1   5798.3   5012.5;
            6364.0       0       3724.3   6649.1   8750.0   4301.2   9250.0   6500.0   7630.8   8450.0   8944.3;
            2651.9    3724.3     0       3201.6   5200.0   1700.0   5700.0   4949.7   3605.6   6150.0   5560.1;
            3162.3    6649.1   3201.6     0       2061.6   4200.0   2500.0   6500.0   1100.0   6100.0   3316.6;
            5000.0    8750.0   5200.0   2061.6     0       6200.0   1118.0   8500.0   3162.3   5830.9   3354.1;
            2702.9    4301.2   1700.0   4200.0   6200.0     0       6708.2   2236.1   3719.5   4450.0   3600.0;
            5025.0    9250.0   5700.0   2500.0   1118.0   6708.2     0       9000.0   3605.6   4031.1   1500.0;
            6324.6    6500.0   4949.7   6500.0   8500.0   2236.1   9000.0     0       5920.4   1950.0   5550.0;
            1581.1    7630.8   3605.6   1100.0   3162.3   3719.5   3605.6   5920.4     0       5250.0   3100.0;
            5798.3    8450.0   6150.0   6100.0   5830.9   4450.0   4031.1   1950.0   5250.0     0       4250.0;
            5012.5    8944.3   5560.1   3316.6   3354.1   3600.0   1500.0   5550.0   3100.0   4250.0     0
        ];
        % 客户点时间窗矩阵（仅用于约束和时间成本计算）
        TW_sec = [
            20*60, 40*60;   % 客户点1
            10*60, 30*60;   % 客户点2
             5*60, 25*60;   % 客户点3
             0*60, 20*60;   % 客户点4
             5*60, 25*60;   % 客户点5
            20*60, 40*60;   % 客户点6
            15*60, 35*60;   % 客户点7
            10*60, 30*60;   % 客户点8
             0*60, 20*60;   % 客户点9
             5*60, 25*60    % 客户点10
             0,100000000000 % 配送中心
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
        Penalty1=0; Penalty2=0; 
        route=PopDec(j,:);
        
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
        
            % 当前剩余载荷（从前往后送）
            segforload=[seg n+1];
            load_kg = zeros(npts-1,1);
            for i = 1:npts-1 %
                remaining = segforload(i:end);  % 还未配送的客户
                load_kg(i) = sum(demand_q(remaining));
            end
            load_all(k) = load_kg(1);% 保存每架无人机初始载荷，以计算约束
        
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

               % ---- 基于“到达点”的时间窗计算时间成本 ----
               % 当前段是 segforlength(i) -> segforlength(i+1)
               % 到达点索引（0 代表仓库，其余为客户 1..n）
               dest_idx = segforlength(i+1);
               if dest_idx == 0
                   % 仓库对应 TW_sec 的最后一行（n+1）
                   tw_row = n + 1;
               else
                   tw_row = dest_idx;
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
        
        PopObj(j,1)=PopObj(j,1)+10*(Penalty1+Penalty2);  % 降低惩罚系数从10^3到10
        PopObj(j,2)=PopObj(j,2)+10*(Penalty1+Penalty2);  % 降低惩罚系数从10^3到10
        end
        
    else
                    j=1;
        Penalty1=0; Penalty2=0; 
        route=PopDec(j,:);
        
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
        
            % 当前剩余载荷（从前往后送）
            segforload=[seg n+1];
            load_kg = zeros(npts-1,1);
            for i = 1:npts-1 %
                remaining = segforload(i:end);  % 还未配送的客户
                load_kg(i) = sum(demand_q(remaining));
            end
            load_all(k) = load_kg(1);% 保存每架无人机初始载荷，以计算约束
        
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

               % ---- 基于“到达点”的时间窗计算时间成本 ----
               % 当前段是 segforlength(i) -> segforlength(i+1)
               % 到达点索引（0 代表仓库，其余为客户 1..n）
               dest_idx = segforlength(i+1);
               if dest_idx == 0
                   % 仓库对应 TW_sec 的最后一行（n+1）
                   tw_row = n + 1;
               else
                   tw_row = dest_idx;
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
            fprintf('解 %d: 惩罚值: Penalty1=%.2f, Penalty2=%.2f\n', j, Penalty1, Penalty2);
        end
        
        PopObj(j,1)=PopObj(j,1)+10*(Penalty1+Penalty2);  % 降低惩罚系数从10^3到10
        PopObj(j,2)=PopObj(j,2)+10*(Penalty1+Penalty2);  % 降低惩罚系数从10^3到10
        
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


