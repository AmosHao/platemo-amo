classdef planner_order_v3 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>

    properties
        % 控制距离矩阵计算方式：
        %   true  → 欧氏距离 + 障碍物惩罚（禁飞区/人流区/噪音区，默认）
        %   false → 纯欧氏距离（忽略所有地图惩罚，方便对比）
        %
        % 切换方式：在 runplatemo_amo.m 调用 platemo() 之前设置全局变量：
        %   global AMOS_USE_MAP_PENALTY;
        %   AMOS_USE_MAP_PENALTY = false;   % 切换为纯欧氏距离
        %   AMOS_USE_MAP_PENALTY = true;    % 恢复含惩罚（或 clear global）
        use_map_penalty = true;
    end

 methods
        %% Default settings of the problem
        function Setting(obj)
            % 检查全局变量，允许在运行脚本里灵活切换距离计算模式
            global AMOS_USE_MAP_PENALTY;
            if ~isempty(AMOS_USE_MAP_PENALTY)
                obj.use_map_penalty = logical(AMOS_USE_MAP_PENALTY);
            end
            obj.n=10;  % 客户点数量
            obj.m=3;   % 无人机数量
            if isempty(obj.M); obj.M = 2; end
            % 编码维度：2*n个点（n个商家+n个客户点）+ (m-1)个0分隔符
            if isempty(obj.D); obj.D = 2*obj.n + (obj.m - 1); end
            obj.lower = zeros(1,obj.D);  % 允许0（分隔符）和1-20（商家和客户点）
            obj.upper = 2*obj.n*ones(1,obj.D);  % 最大编号为20（商家编号1-10，客户点编号11-20）
            % 注意：不能使用encoding=5（permutation），因为序列中包含重复的0分隔符
            % 使用encoding=2（integer），然后通过自定义操作符处理
            obj.encoding = 2*ones(1,obj.D);  % 改为integer编码
            % 增加种群大小以提高搜索多样性和避免局部最优
            % 排列组合空间分析：
            % - 20个点（商家1-10，客户点11-20）需要访问
            % - 需要分成3段（3架无人机），分段方式：C(19,2) = 171种
            % - 考虑顺序约束后，排列组合空间仍然很大（远大于20!）
            % 因此需要足够大的种群来充分探索解空间
            if isempty(obj.N) || obj.N < 300
                obj.N = 300;  % 从200增加到300，进一步提高多样性
            end
        end
        %% Generate initial solutions
        function Population = Initialization(obj,N)
            if nargin < 2
                N = obj.N;
            end
            n = obj.n;  % 客户点数量（商家数量）
            m = obj.m;  % 无人机数量
            D = obj.D;  % 决策变量维度：2*n + (m-1)
            
            % 初始化 PopDec 矩阵
            PopDec = zeros(N, D);
            
            % 遍历每一行生成初始解
            for i = 1:N
                % 策略：生成满足顺序约束的初始解
                % 对于每个商家-客户点对，确保商家在客户点之前
                
                % 生成商家-客户点对（每个对包含一个商家和对应客户点）
                pairs = cell(n, 1);
                for p = 1:n
                    % 每个对：商家p和客户点p+n
                    pairs{p} = [p, p+n];
                end
                
                % 随机打乱对的顺序
                pair_order = randperm(n);
                
                % 构建满足顺序约束的序列
                % 确保所有初始解都满足顺序约束：商家必须在对应客户点之前
                % 多样性通过以下方式实现：
                % 1. 不同的商家-客户点对排列顺序（pair_order）
                % 2. 不同的0分隔符位置
                % 3. 不同的段内顺序（某些对可以紧邻，某些对可以分开）
                sequence = [];
                
                for p_idx = 1:length(pair_order)
                    pair = pairs{pair_order(p_idx)};
                    % 始终满足顺序约束：商家在前，客户在后
                    sequence = [sequence, pair(1), pair(2)];  % 商家，客户
                end
                
                % 随机选择 m-1 个位置插入0分隔符
                % 关键：由于商家-客户点是成对出现的（每个对2个元素），
                % 0分隔符必须插入在偶数位置之后（客户点之后），不能插入在奇数位置（商家之后），
                % 否则会拆散商家-客户点对
                if m > 1 && length(sequence) > 1
                    % sequence长度是2*n（n个商家-客户点对）
                    % 可插入位置：2, 4, 6, ..., 2*n-2（偶数位置，客户点之后）
                    % 这样不会拆散商家-客户点对
                    even_positions = 2:2:(2*n-2);  % 所有可用的偶数位置
                    
                    if length(even_positions) >= m - 1
                        % 从偶数位置中随机选择 m-1 个
                        selected_indices = randperm(length(even_positions), m - 1);
                        positions = even_positions(selected_indices);
                        positions = sort(positions);  % 排序
                        
                        % 在选定的位置后插入0
                        result = [];
                        prev_pos = 0;
                        for j = 1:length(positions)
                            pos = positions(j);
                            result = [result, sequence(prev_pos+1:pos), 0];
                            prev_pos = pos;
                        end
                        result = [result, sequence(prev_pos+1:end)];
                    else
                        % 如果偶数位置不足，在末尾补充0
                        result = sequence;
                        result = [result, zeros(1, m - 1)];
                    end
                else
                    result = sequence;
                end
                
                % 检查长度是否正确（用于调试，不执行修复）
                if length(result) ~= D
                    % 打印详细错误信息
                    fprintf('\n========== 初始化解长度错误 [解 %d] ==========\n', i);
                    fprintf('期望长度 D = %d (2*n + (m-1) = 2*%d + (%d-1) = %d)\n', D, n, m, D);
                    fprintf('实际长度 = %d\n', length(result));
                    fprintf('长度差异 = %d\n', length(result) - D);
                    fprintf('0分隔符数量 = %d (期望: %d)\n', sum(result == 0), m - 1);
                    fprintf('非0值数量 = %d (期望: %d)\n', sum(result ~= 0), 2*n);
                    fprintf('sequence原始长度 = %d\n', length(sequence));
                    fprintf('result前20个元素: [%s]\n', num2str(result(1:min(20, length(result)))));
                    if length(result) > 20
                        fprintf('result后20个元素: [%s]\n', num2str(result(end-min(19, length(result)-20):end)));
                    end
                    fprintf('==========================================\n\n');
                    % 抛出错误，停止执行
                    error('初始化解长度不匹配：期望长度=%d，实际长度=%d', D, length(result));
                end
                
                % 最终验证：检查0分隔符数量是否正确（用于调试，不执行修复）
                zeroCount = sum(result == 0);
                if zeroCount ~= m - 1
                    % 打印详细错误信息
                    fprintf('\n========== 初始化解0分隔符数量错误 [解 %d] ==========\n', i);
                    fprintf('期望0分隔符数量 = %d (m-1 = %d-1 = %d)\n', m - 1, m, m - 1);
                    fprintf('实际0分隔符数量 = %d\n', zeroCount);
                    fprintf('数量差异 = %d\n', zeroCount - (m - 1));
                    fprintf('result长度 = %d (期望: %d)\n', length(result), D);
                    zeroIdx = find(result == 0);
                    if ~isempty(zeroIdx)
                        fprintf('0分隔符位置: [%s]\n', num2str(zeroIdx));
                    else
                        fprintf('0分隔符位置: 无\n');
                    end
                    fprintf('result前20个元素: [%s]\n', num2str(result(1:min(20, length(result)))));
                    if length(result) > 20
                        fprintf('result后20个元素: [%s]\n', num2str(result(end-min(19, length(result)-20):end)));
                    end
                    fprintf('==========================================\n\n');
                    % 抛出错误，停止执行
                    error('初始化解0分隔符数量错误：期望=%d，实际=%d', m - 1, zeroCount);
                end
                
                PopDec(i, :) = result;
            end
            
            % 评估初始种群，用掉N次
            Population = obj.Evaluation(PopDec); 
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
        % 空气密度和旋翼桨盘面积（用于计算诱导速度 v_0）
        rho_air = 1.225;  % 空气密度 (kg/m³)
        A_rotor = 0.503;  % 旋翼桨盘面积 (m²)
        % 论文参数（表2-1）
        Omega = 300;      % 旋翼角速度 (rad/s)
        R_rotor = 0.4;    % 旋翼半径 (m)
        delta_drag = 0.012;  % 平均阻力系数
        kappa = 1.1;      % 诱导功率修正系数
        S_FP = 0.0151;   % 等效机身面积 (m²)
        % 计算型阻功率常数 P_0 = (δ/8)ρsAΩ³R³
        % 假设桨盘实度 s = bcR/(πR²)，对于4旋翼，b=4，c=0.0157
        b_blades = 4;    % 旋翼数量
        c_chord = 0.0157; % 翼型弦长 (m)
        s_solidity = b_blades * c_chord * R_rotor / (pi * R_rotor^2);  % 桨盘实度
        P_0 = (delta_drag / 8) * rho_air * s_solidity * A_rotor * Omega^3 * R_rotor^3;  % 型阻功率常数
        % 引入数据文件（包含坐标、障碍物、惩罚系数等）
        % 数据文件路径：从当前文件位置到 forOrderNew26/Order_Map/order_data.m
        current_file_dir = fileparts(mfilename('fullpath'));
        % 从 PlatEMO/Problems/Multi-objective optimization/Amos/ 到 PlatEMO/ 目录
        platemo_root = fileparts(fileparts(fileparts(current_file_dir)));
        data_file_path = fullfile(platemo_root, 'forOrderNew26', 'Order_Map', 'order_data.m');
        
        if exist(data_file_path, 'file')
            run(data_file_path);
        else
            % 如果找不到数据文件，使用默认数据（向后兼容）
            warning('找不到数据文件 %s，使用默认数据', data_file_path);
            % // dotss = [
            % //     5000, 5000, 800;   % 配送中心（在正中间）
            % //     8156, 1270, 800;   % 商家1
            % //     9058, 2785, 800;   % 商家2
            % //     8500, 1576, 800;   % 商家3
            % //      975, 1419, 800;   % 商家4
            % //     2785, 5469, 800;   % 商家5
            % //     9575, 9649, 800;   % 商家6
            % //     4854, 8003, 800;   % 商家7
            % //     7922, 9595, 800;   % 商家8
            % //      357, 6557, 800;   % 商家9
            % //     8491, 9340, 800;   % 商家10
            % //      500,  500, 800;   % 客户点11
            % //     1500, 2800, 800;   % 客户点12
            % //     7000, 3000, 800;   % 客户点13
            % //     9000, 5500, 800;   % 客户点14
            % //     1500, 4500, 800;   % 客户点15
            % //     9500, 5500, 800;   % 客户点16
            % //      500, 7000, 800;   % 客户点17
            % //     6000, 6100, 800;   % 客户点18
            % //     1800, 8950, 800;   % 客户点19
            % //     6050, 9200, 800    % 客户点20
            % // ];
            % // forbidden_zones = [
            % //     2000, 2000, 4000, 4000;
            % //     6000, 1000, 8000, 3000;
            % // ];
            % // crowded_zones = [
            % //     6500, 6500, 500;
            % //     3000, 7000, 800;
            % // ];
            % // noise_zones = [
            % //     4500, 3500, 5500, 4500;
            % //     8500, 3500, 9500, 4500;
            % // ];
            % // penalty_forbidden = 2000;
            % // penalty_crowded = 500;
            % // penalty_noise = 1000;
        end
        
        % 计算距离矩阵（21×21）
        % 行/列索引：1-10=商家1-10, 11-20=客户点11-20, 21=配送中心
        if obj.use_map_penalty
            % 含障碍物惩罚（禁飞区/人流区/噪音区）
            d_matrix = calculate_distance_matrix_with_obstacles(...
                dotss, forbidden_zones, crowded_zones, noise_zones, ...
                penalty_forbidden, penalty_crowded, penalty_noise);
        else
            % 纯欧氏距离（传入空障碍物，函数会跳过所有惩罚）
            d_matrix = calculate_distance_matrix_with_obstacles(...
                dotss, [], [], [], 0, 0, 0);
        end
        % 时间窗矩阵（21行：商家1-10，客户点11-20，配送中心）
        % 行索引：1-10=商家1-10, 11-20=客户点11-20, 21=配送中心
        TW_sec = [
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
             5*60, 25*60;     % 客户点20
            0, 100000000000   % 配送中心（第21行）
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
        % 修正后的能量公式 + 10km×10km地图：
        %   均衡解（每架3-4单,7跳,均距4km）：~540kJ/架 → 900kJ 时可行
        %   偏斜解（某架5单,11跳）：~850kJ → 900kJ 时仍可行（约束有余量）
        %   极端解（某架7单以上）：>1000kJ → 仍违约（约束起惩罚作用）
        % 结论：10km地图下900kJ是使约束有意义的合理值（约250Wh电池）
        maxEC=800000;%电池容量（J），800kJ

    if size(PopDec,1)>0
        % 调试统计
    debug_stats = struct();
    debug_stats.invalid_reasons = struct();
    debug_stats.invalid_reasons.inf_nan = 0;      % Inf或NaN
    debug_stats.invalid_reasons.negative = 0;     % 负数
    debug_stats.invalid_reasons.missing_points = 0;  % 缺少必需的点
    debug_stats.invalid_reasons.wrong_zeros = 0;  % 0分隔符数量错误
    debug_stats.invalid_reasons.empty_segments = 0;  % 空段
    debug_stats.invalid_reasons.energy_error = 0;  % 能量计算错误
    debug_stats.invalid_reasons.time_error = 0;   % 时间计算错误
    debug_stats.total_invalid = 0;
    debug_stats.invalid_examples = {};  % 保存前几个无效解的示例
    
    for j = 1:N 
        Penalty1=0; Penalty2=0; Penalty3=0;  % Penalty3用于顺序约束
        route=PopDec(j,:);
        
        % 调试：检查解的结构
        nonZeroVals = route(route ~= 0);
        zeroCount = sum(route == 0);
        allRequiredPoints = 1:(2*n);
        missingPoints = setdiff(allRequiredPoints, unique(nonZeroVals));
        
        % 检查解的结构完整性
        if ~isempty(missingPoints)
            debug_stats.invalid_reasons.missing_points = debug_stats.invalid_reasons.missing_points + 1;
            if length(debug_stats.invalid_examples) < 3
                debug_stats.invalid_examples{end+1} = struct('type', 'missing_points', ...
                    'route', route, 'missing', missingPoints, 'zeroCount', zeroCount);
            end
        end
        
        if zeroCount ~= (m-1)
            debug_stats.invalid_reasons.wrong_zeros = debug_stats.invalid_reasons.wrong_zeros + 1;
            if length(debug_stats.invalid_examples) < 3
                debug_stats.invalid_examples{end+1} = struct('type', 'wrong_zeros', ...
                    'route', route, 'zeroCount', zeroCount, 'expected', m-1);
            end
        end
        
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
        % 检查0分隔符数量是否正确（用于调试，不执行修复）
        if length(idx0) ~= m - 1
            fprintf('计算目标函数时检测到0的数量不对 [解 %d]: 期望=%d, 实际=%d\n', j, m - 1, length(idx0));
        end
        % 使用实际的0分隔符数量进行分段（即使数量不对也继续计算）
        idx0   = [0 idx0 numel(route)+1];          % 方便循环
        routes = cell(m,1);
        for k = 1:m
            if k <= length(idx0) - 1
                routes{k} = route(idx0(k)+1 : idx0(k+1)-1);
            else
                routes{k} = [];  % 如果段数不足，返回空
            end
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
                debug_stats.invalid_reasons.empty_segments = debug_stats.invalid_reasons.empty_segments + 1;
                if length(debug_stats.invalid_examples) < 3
                    debug_stats.invalid_examples{end+1} = struct('type', 'empty_segments', ...
                        'route', route, 'emptySeg', k);
                end
                continue;
            end
        
            % 路径点：depot -> seg(1) -> seg(2) -> ... -> seg(end) -> depot
            npts = length(seg) + 2; % 段内点数+2（depot起点和终点）
        
            % 当前剩余载荷计算（考虑商家取餐和客户点配送）
            % 关键：load_kg(i) 应该表示从 segforlength(i) 起飞时的载重
            % segforlength = [0, seg, 0] 表示：depot -> seg(1) -> seg(2) -> ... -> seg(end) -> depot
            segforlength = [0,seg,0];
            % 计算每个点的载重状态（到达该点后的载重）
            load_at_point = zeros(length(segforlength), 1);
            load_at_point(1) = 0;  % depot出发时载重为0
            
            for i = 2:length(segforlength)-1
                point_id = segforlength(i);
                if point_id >= 1 && point_id <= n
                    % 商家（编号1-10）：取餐，载重增加
                    customer_idx = point_id;  % 商家i对应客户点i+10，但demand_q索引是1-10
                    load_at_point(i) = load_at_point(i-1) + demand_q(customer_idx);
                elseif point_id > n && point_id <= 2*n
                    % 客户点（编号11-20）：配送，载重减少
                    customer_idx = point_id - n;  % 客户点编号转换为demand_q索引（1-10）
                    load_at_point(i) = load_at_point(i-1) - demand_q(customer_idx);
                else
                    load_at_point(i) = load_at_point(i-1);  % 其他情况（不应该发生）
                end
            end
            load_at_point(end) = 0;  % 回到depot时载重为0
            
            % load_kg(i) 表示从 segforlength(i) 起飞时的载重（即到达 segforlength(i) 后的载重）
            load_kg = load_at_point(1:end-1);  % 去掉最后一个点（终点depot）
            load_all(k) = max(load_at_point);  % 保存每架无人机最大载重，以计算约束
        
               T_k_all = 0; 
               E_k     = 0;
               % 累积时间（到达每个点的时间）
               time_cum = 0;
            
            for i = 1:npts-1 
                % 此处引入航迹规划算法，输入dot索引和坐标, 输出路程。需要缓存判断。
                % [length]=runPlatemo_forOrder_PRMO(problem,algorithm,N,maxFE,s,M,dots,dot1,dot2,whichObj);
                % 索引映射：0=配送中心(21), 1-10=商家(1-10), 11-20=客户点(11-20)
                if segforlength(i) == 0
                    dot1 = 21;  % 配送中心
                else
                    dot1 = segforlength(i);  % 商家或客户点
                end
                if segforlength(i+1) == 0
                    dot2 = 21;  % 配送中心
                else
                    dot2 = segforlength(i+1);  % 商家或客户点
                end 
                % dots=[dotss(dot1,:),dotss(dot2,:),];
                % if d_matrix(dot1,dot2)==0
                % [f_length]=runPlatemo_forOrder_PRMO(@planner_simple_maxhjd_newobj,@PRMO,100,300,20,3,dots,dot1,dot2,1);
                % d_matrix(dot1,dot2)=f_length;
                % d_matrix(dot2,dot1)=f_length;
                % d_hor=f_length;
                % else
                d_hor=d_matrix(dot1,dot2);
                % end
                
                % 调试：检查距离矩阵索引
                if dot1 < 1 || dot1 > size(d_matrix,1) || dot2 < 1 || dot2 > size(d_matrix,2)
                    fprintf('错误 [解 %d, 无人机 %d, 段 %d]: 距离矩阵索引越界 - dot1=%d, dot2=%d, 矩阵大小=[%d,%d]\n', ...
                        j, k, i, dot1, dot2, size(d_matrix,1), size(d_matrix,2));
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                    if length(debug_stats.invalid_examples) < 5
                        debug_stats.invalid_examples{end+1} = struct('type', 'invalid_dot_index', ...
                            'route', route, 'dot1', dot1, 'dot2', dot2, 'matrixSize', size(d_matrix), ...
                            'seg', seg, 'i', i, 'k', k);
                    end
                end
                
                % 调试：检查距离值
                if ~isfinite(d_hor) || d_hor < 0
                    fprintf('错误 [解 %d, 无人机 %d, 段 %d]: 距离值无效 - d_hor=%f\n', j, k, i, d_hor);
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                end

                % 上升（从所有点起飞时都需要上升，因为到达时都会降落）
                t_up = h_up_down / v_u;
                
                % 调试：检查时间计算
                if ~isfinite(t_up) || t_up < 0 || v_u <= 0
                    debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
                    t_up = max(0.1, t_up);  % 使用最小值
                end
                % 此处写入上升功率计算公式
                % load_kg(i) 现在正确表示从 segforlength(i) 起飞时的载重
                UAV_w = (UAV_m+load_kg(i))*g; %无人机重量（N）
                % 诱导速度 v_0 = sqrt((m+q)g / (2*rho*A))
                % 其中：m+q 是总质量（kg），g 是重力加速度，rho 是空气密度，A 是旋翼桨盘面积
                UAV_mass = UAV_m + load_kg(i);  % 总质量（kg）
                v_0 = sqrt((UAV_mass * g) / (2 * rho_air * A_rotor));
                sqrt_term_up = sqrt((v_u/2)^2 + v_0^2);  % 量纲修正：括号内应为 v_0^2（m²/s²），而非 v_0（m/s）
                P_up_i = P_0 + UAV_w*(v_u/2 + sqrt_term_up);  % P_0 代替硬编码常数，与动态计算一致
                E_up = P_up_i * t_up;
        
                % 水平飞行
                t_hor = d_hor / v_h;%此处d_hor为航迹规划算法输出的值或矩阵内存调用
                % 水平功率计算公式（根据论文公式2-38）
                % P_h = P_0(1 + 3μ²) + P_i * sqrt_term + (1/2)ρS_FP*V_h³
                % 其中：μ = V_h/(ΩR)，P_i = κ*T²/sqrt(2ρA)
                
                % 1. 型阻功率：P_pro,h = P_0(1 + 3μ²)
                mu = v_h / (Omega * R_rotor);  % 前进比
                P_pro_h = P_0 * (1 + 3 * mu^2);
                
                % 2. 诱导功率：P_i,h = κ*T²/sqrt(2ρA) * sqrt_term
                % sqrt_term = (sqrt(1 + V_h⁴/(4v_0⁴)) - V_h²/(2v_0²))^(1/2)
                % 正确公式（Zeng et al.）：sqrt(1 + V^4/(4*v_0^4)) - V^2/(2*v_0^2)
                % 注意：分子内的 1 独立于 /(4*v_0^4)，不能写成 (1+v_h^4)/(4*v_0^4)（否则两项几乎相消，P_i_h≈0）
                sqrt_inner = sqrt(1 + v_h^4/(4*v_0^4)) - v_h^2/(2*v_0^2);
                sqrt_term_hor = sqrt(max(0, sqrt_inner));  % max(0,·) 防止数值误差产生负数
                % P_i_h = κ·T^(3/2)/sqrt(2ρA)·sqrt_term = κ·T·v_0·sqrt_term（T^1.5 而非 T^2）
                P_i_h = kappa * UAV_w * v_0 * sqrt_term_hor;
                
                % 3. 废阻功率：P_par = (1/2)ρS_FP*V_h³
                P_par = 0.5 * rho_air * S_FP * v_h^3;
                
                % 总水平功率
                P_hor_i = P_pro_h + P_i_h + P_par;
                E_hor = P_hor_i * t_hor;
        
                % 下降（到达所有点时都需要下降，包括depot、商家和客户点）
                t_down = h_up_down / v_d;
                %此处写入下降功率计算公式
                v_d_0=v_d/v_0;
                P_down_i = 79.85625 + UAV_w*v_0*(0.974 -1.125*v_d_0 -1.372*v_d_0^2 -1.718*v_d_0^3 -0.655*v_d_0^4);
                % 注意：在风车状态（autorotation）下，P_down_i 可能为负（表示能量回收或几乎不消耗能量）
                % 由于大多数无人机没有能量回收系统，负功率时能量消耗设为0
                E_down = max(0, P_down_i * t_down);
                
                % 检测复数或异常功率（标记为无效，不使用real()掩盖）
                % 注意：修正后 sqrt_inner 应该总是为正
                % 注意：P_down_i < 0 在风车状态下是合理的，不再标记为异常
                has_complex_or_negative = false;
                if ~isreal(sqrt_term_up) || (sqrt_inner < 0 || ~isreal(sqrt_term_hor)) || ...
                   ~isreal(P_up_i) || ~isreal(E_up) || ...
                   ~isreal(P_hor_i) || ~isreal(E_hor) || ...
                   ~isreal(P_down_i) || ~isreal(E_down)
                    has_complex_or_negative = true;
                    % 计算 sqrt_inner 的两个组成部分用于调试
                    term1 = sqrt((1+v_h^4)/(4*v_0^4));
                    term2 = v_h^2/(2*v_0^2);  % 修正：应该是除以 2*v_0^2
                    fprintf('警告 [解 %d, 无人机 %d, 段 %d]: 功率计算异常（将标记为无效解）\n', j, k, i);
                    fprintf('  v_0=%.4f, UAV_mass=%.2f, load_kg=%.2f, UAV_w=%.2f\n', v_0, UAV_mass, load_kg(i), UAV_w);
                    fprintf('  sqrt_inner=%.4f (term1=%.4f, term2=%.4f, term1-term2=%.4f)\n', sqrt_inner, term1, term2, term1-term2);
                    fprintf('  sqrt_term_hor=%s, P_hor_i=%s\n', num2str(sqrt_term_hor), num2str(P_hor_i));
                    fprintf('  P_down_i=%.4f, v_d_0=%.4f, v_d=%.1f\n', P_down_i, v_d_0, v_d);
                    % 将能量和时间设置为NaN，后续会被标记为无效解
                    E_up = NaN;
                    E_hor = NaN;
                    E_down = NaN;
                end

               % 计算去往每个点路上所花费的时间（单段）
               T_k_single = t_up + t_hor + t_down;
               % 累积到达时间（从出发到当前“到达点”的总时间）
               time_cum = time_cum + T_k_single;

               % ---- 基于"到达点"的时间窗计算时间成本 ----
               % 当前段是 segforlength(i) -> segforlength(i+1)
               % 到达点索引：0=配送中心(21), 1-10=商家(1-10), 11-20=客户点(11-20)
               dest_idx = segforlength(i+1);
               if dest_idx == 0
                   % 配送中心对应 TW_sec 的第21行
                   tw_row = 21;
               elseif dest_idx >= 1 && dest_idx <= n
                   % 商家（1-10）对应 TW_sec 的第1-10行
                   tw_row = dest_idx;
               elseif dest_idx > n && dest_idx <= 2*n
                   % 客户点（11-20）对应 TW_sec 的第11-20行
                   tw_row = dest_idx;  % 客户点11对应第11行，客户点20对应第20行
               else
                   tw_row = 21;  % 默认使用配送中心时间窗
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
                % 调试：检查累加前的值
                if ~isfinite(T_k) || T_k < 0
                    debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
                end
                if ~isfinite(E_up) || ~isfinite(E_hor) || ~isfinite(E_down)
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                end
                
                T_k_all = T_k_all + T_k;
                E_k = E_k + E_up + E_hor + E_down;
                % 确保累加后的值也是实数
                T_k_all = real(T_k_all);
                E_k = real(E_k);
                
                % 调试：检查累加后的值
                if ~isfinite(T_k_all) || ~isfinite(E_k)
                    debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                end
            end
        
            T_all(k) = real(T_k_all);%每架飞机的
            E_all(k) = real(E_k);
            
            % 调试：检查最终的能量和时间值
            if ~isfinite(T_all(k)) || T_all(k) < 0
                debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
            end
            if ~isfinite(E_all(k)) || E_all(k) < 0
                debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
            end
        end
        
        % 目标一：总能耗（移除区分度项，直接使用总和）
        % 确保 E_all 和 T_all 都是实数
        E_all = real(E_all);
        T_all = real(T_all);
        
        baseE_sum = real(sum(E_all));
        PopObj(j,1) = baseE_sum;
        PopObj(j,1) = real(PopObj(j,1));
        
        % 目标二：总时间（移除区分度项，直接使用总和）
        baseT_sum = real(sum(T_all));
        PopObj(j,2) = baseT_sum;
        PopObj(j,2) = real(PopObj(j,2)); 
         
        % 约束--载重
        for kk = 1:m
            load_k = load_all(kk);
            if load_k > maxload
                Penalty1=Penalty1+(load_k-maxload);
            end
        end
        
        % 约束--电量（直接用焦耳，与 Obj1 同单位；不再除以 1000）
        for kk = 1:m
            power_k = E_all(kk);
            if power_k > maxEC
                Penalty2 = Penalty2 + (power_k - maxEC);
            end
        end
        
        % ---- 惩罚项加到目标函数 ----
        % 设计原则：
        %   Obj1（总能耗，J）：违约量本身就是 J，直接加，无需缩放
        %   Obj2（时间成本，s）：用固定物理换算 1J ≈ 1/150s（约150W平均飞行功率）
        %       使时间目标中的能量惩罚与时间值量级相当（100kJ→667s，典型T约1000-3000s）
        %   两者使用相同的违约量 → 无论个体能量高低，惩罚比例一致，不扭曲Pareto前沿
        E_to_T_scale = 1/150;  % 固定物理换算比（s/J），基于约150W平均水平巡航功率
        totalViol = Penalty1 + Penalty2;  % 总违约量（J当量，Penalty1已乘loadScale=1e4）
        
        PopObj(j,1) = real(PopObj(j,1) + totalViol);
        PopObj(j,2) = real(PopObj(j,2) + totalViol * E_to_T_scale + Penalty3);
        
        % 确保目标值为有限值（记录无效原因）
        is_invalid = false;
        invalid_reason = '';
        
        if ~isfinite(PopObj(j,1)) || isnan(PopObj(j,1))
            is_invalid = true;
            invalid_reason = 'f1_Inf_NaN';
            debug_stats.invalid_reasons.inf_nan = debug_stats.invalid_reasons.inf_nan + 1;
        elseif PopObj(j,1) < 0
            is_invalid = true;
            invalid_reason = 'f1_negative';
            debug_stats.invalid_reasons.negative = debug_stats.invalid_reasons.negative + 1;
        end
        
        if ~isfinite(PopObj(j,2)) || isnan(PopObj(j,2))
            is_invalid = true;
            if isempty(invalid_reason)
                invalid_reason = 'f2_Inf_NaN';
            else
                invalid_reason = [invalid_reason, '_f2_Inf_NaN'];
            end
            debug_stats.invalid_reasons.inf_nan = debug_stats.invalid_reasons.inf_nan + 1;
        elseif PopObj(j,2) < 0
            is_invalid = true;
            if isempty(invalid_reason)
                invalid_reason = 'f2_negative';
            else
                invalid_reason = [invalid_reason, '_f2_negative'];
            end
            debug_stats.invalid_reasons.negative = debug_stats.invalid_reasons.negative + 1;
        end
        
        if is_invalid
            PopObj(j,1) = 1e10;  % 设置为一个大的有限值
            PopObj(j,2) = 1e10;  % 设置为一个大的有限值
            debug_stats.total_invalid = debug_stats.total_invalid + 1;
            
            % 保存无效解示例
            if length(debug_stats.invalid_examples) < 5
                debug_stats.invalid_examples{end+1} = struct('type', invalid_reason, ...
                    'route', route, 'E_all', E_all, 'T_all', T_all, ...
                    'PopObj_before', [baseE_sum, baseT_sum], 'Penalty1', Penalty1, ...
                    'Penalty2', Penalty2, 'Penalty3', Penalty3, 'PopObj_after', PopObj(j,:));
            end
        end
        end
        
        % 输出调试统计信息（每处理完所有解后输出一次）
        % if j == N
        %     fprintf('\n========== 调试统计 [批次完成] ==========\n');
        %     fprintf('总解数: %d\n', N);
        %     fprintf('无效解总数: %d (%.1f%%)\n', debug_stats.total_invalid, 100*debug_stats.total_invalid/N);
        %     fprintf('\n无效原因统计:\n');
        %     fprintf('  - Inf/NaN: %d\n', debug_stats.invalid_reasons.inf_nan);
        %     fprintf('  - 负数: %d\n', debug_stats.invalid_reasons.negative);
        %     fprintf('  - 缺少必需的点: %d\n', debug_stats.invalid_reasons.missing_points);
        %     fprintf('  - 0分隔符数量错误: %d\n', debug_stats.invalid_reasons.wrong_zeros);
        %     fprintf('  - 空段: %d\n', debug_stats.invalid_reasons.empty_segments);
        %     fprintf('  - 能量计算错误: %d\n', debug_stats.invalid_reasons.energy_error);
        %     fprintf('  - 时间计算错误: %d\n', debug_stats.invalid_reasons.time_error);
        % 
        %     if ~isempty(debug_stats.invalid_examples)
        %         fprintf('\n无效解示例（前%d个）:\n', min(3, length(debug_stats.invalid_examples)));
        %         for ex_idx = 1:min(3, length(debug_stats.invalid_examples))
        %             ex = debug_stats.invalid_examples{ex_idx};
        %             fprintf('  示例%d - 类型: %s\n', ex_idx, ex.type);
        %             fprintf('    路径长度: %d, 路径前15个: [%s]\n', length(ex.route), ...
        %                 num2str(ex.route(1:min(15, length(ex.route)))));
        %             if isfield(ex, 'missing')
        %                 fprintf('    缺失的点: [%s]\n', num2str(ex.missing));
        %             end
        %             if isfield(ex, 'zeroCount')
        %                 if isfield(ex, 'expected')
        %                     fprintf('    0分隔符数量: %d (期望: %d)\n', ex.zeroCount, ex.expected);
        %                 else
        %                     fprintf('    0分隔符数量: %d (期望: %d)\n', ex.zeroCount, m-1);
        %                 end
        %             end
        %             if isfield(ex, 'E_all')
        %                 fprintf('    E_all: [%s]\n', num2str(ex.E_all));
        %             end
        %             if isfield(ex, 'T_all')
        %                 fprintf('    T_all: [%s]\n', num2str(ex.T_all));
        %             end
        %             if isfield(ex, 'PopObj_before')
        %                 fprintf('    原始目标值: f1=%.2f, f2=%.2f\n', ex.PopObj_before(1), ex.PopObj_before(2));
        %             end
        %             if isfield(ex, 'Penalty1')
        %                 fprintf('    惩罚值: P1=%.2f, P2=%.2f, P3=%.2f\n', ex.Penalty1, ex.Penalty2, ex.Penalty3);
        %             end
        %         end
        %     end
        %     fprintf('==========================================\n\n');
        % end
        
    else

        fprintf('PopDec is empty\n');
        %             j=1;
        % Penalty1=0; Penalty2=0; Penalty3=0;  % Penalty3用于顺序约束
        % route=PopDec(j,:);
        % 
        % % ---------- 0. 检查顺序约束：商家必须在对应客户点之前访问 ----------
        % % 商家编号：1-10，对应客户点编号：11-20
        % for merchant_id = 1:n
        %     customer_id = merchant_id + n;  % 客户点编号
        %     merchant_pos = find(route == merchant_id);
        %     customer_pos = find(route == customer_id);
        % 
        %     if ~isempty(customer_pos)
        %         if isempty(merchant_pos)
        %             Penalty3 = Penalty3 + 1000;
        %         else
        %             idx0 = find(route == 0);
        %             idx0 = [0 idx0 numel(route)+1];
        % 
        %             customer_seg = 0;
        %             merchant_seg = 0;
        %             for seg_idx = 1:length(idx0)-1
        %                 seg_range = (idx0(seg_idx)+1):(idx0(seg_idx+1)-1);
        %                 if any(route(seg_range) == customer_id)
        %                     customer_seg = seg_idx;
        %                 end
        %                 if any(route(seg_range) == merchant_id)
        %                     merchant_seg = seg_idx;
        %                 end
        %             end
        % 
        %             if merchant_seg == customer_seg && merchant_seg > 0
        %                 seg_range = (idx0(merchant_seg)+1):(idx0(merchant_seg+1)-1);
        %                 merchant_pos_in_seg = find(route(seg_range) == merchant_id, 1);
        %                 customer_pos_in_seg = find(route(seg_range) == customer_id, 1);
        %                 if merchant_pos_in_seg >= customer_pos_in_seg
        %                     Penalty3 = Penalty3 + 1000;
        %                 end
        %             elseif merchant_seg > customer_seg && customer_seg > 0
        %                 Penalty3 = Penalty3 + 1000;
        %             end
        %         end
        %     end
        % end
        % 
        % % ---------- 1. 根据 0 切分m段 ----------
        % idx0   = find(route == 0);
        % % 确保有 m-1 个0分隔符
        % if length(idx0) < m - 1
        %     % 如果0的数量不足，在末尾补充0
        %     numZerosNeeded = (m - 1) - length(idx0);
        %     % 找到非0的位置，在适当位置插入0
        %     nonZeroIdx = find(route ~= 0);
        %     if ~isempty(nonZeroIdx)
        %         % 在非0元素之间插入0
        %         for z = 1:numZerosNeeded
        %             if length(nonZeroIdx) > 1
        %                 insertPos = nonZeroIdx(randi(length(nonZeroIdx)-1)) + 1;
        %                 route = [route(1:insertPos-1), 0, route(insertPos:end)];
        %                 nonZeroIdx = find(route ~= 0);
        %             else
        %                 route = [route, 0];
        %             end
        %         end
        %     else
        %         % 如果全是0，补充到m-1个
        %         route = [route, zeros(1, numZerosNeeded)];
        %     end
        %     idx0 = find(route == 0);
        % elseif length(idx0) > m - 1
        %     % 如果0的数量过多，只保留前m-1个
        %     idx0 = idx0(1:m-1);
        % end
        % idx0   = [0 idx0 numel(route)+1];          % 方便循环
        % routes = cell(m,1);
        % for k = 1:m
        %     if k <= length(idx0) - 1
        %         routes{k} = route(idx0(k)+1 : idx0(k+1)-1);
        %     else
        %         routes{k} = [];  % 如果段数不足，返回空
        %     end
        % end
        % 
        % % 路径按 0 分段后直接用于后续计算
        % % --- 2. 初始化输出 ---
        % T_all = zeros(1,m);  % 每架无人机总时间
        % E_all = zeros(1,m);  % 每架无人机总能耗
        % load_all = zeros(1,m); % 每架无人机初始载重
        % 
        % % --- 3. 每架无人机逐点计算 ---
        % for k = 1:m
        %     seg = routes{k}; % 1 2 3
        %     if isempty(seg)
        %         continue;
        %     end
        % 
        %     % 路径点：depot -> seg(1) -> seg(2) -> ... -> seg(end) -> depot
        %     % pts = [depot; customerXY(seg,:); depot];
        %     npts = length(seg) + 2; %3+2=5
        % 
        %     % 当前剩余载荷计算（考虑商家取餐和客户点配送）
        %     % 商家（1-10）：取餐，载重增加
        %     % 客户点（11-20）：配送，载重减少
        %     segforload = [seg, 2*n+1];  % 添加配送中心作为终点（索引为2*n+1=21）
        %     load_kg = zeros(npts-1, 1);
        %     current_load = 0;  % 从配送中心出发时载重为0
        % 
        %     for i = 1:npts-1
        %         point_id = segforload(i);
        %         if point_id >= 1 && point_id <= n
        %             % 商家（编号1-10）：取餐，载重增加
        %             customer_idx = point_id;  % 商家i对应客户点i+10，但demand_q索引是1-10
        %             current_load = current_load + demand_q(customer_idx);
        %         elseif point_id > n && point_id <= 2*n
        %             % 客户点（编号11-20）：配送，载重减少
        %             customer_idx = point_id - n;  % 客户点编号转换为demand_q索引（1-10）
        %             current_load = current_load - demand_q(customer_idx);
        %         end
        %         load_kg(i) = current_load;
        %     end
        %     load_all(k) = max(load_kg);  % 保存每架无人机最大载重，以计算约束
        % 
        %        T_k_all = 0; 
        %        E_k     = 0;
        %        segforlength = [0,seg,0];
        %        % 累积时间（到达每个点的时间）
        %        time_cum = 0;
        % 
        %     for i = 1:npts-1 
        %         % 此处引入航迹规划算法，输入dot索引和坐标, 输出路程。需要缓存判断。
        %         % [length]=runPlatemo_forOrder_PRMO(problem,algorithm,N,maxFE,s,M,dots,dot1,dot2,whichObj);
        %         % 索引映射：0=配送中心(21), 1-10=商家(1-10), 11-20=客户点(11-20)
        %         if segforlength(i) == 0
        %             dot1 = 21;  % 配送中心
        %         else
        %             dot1 = segforlength(i);  % 商家或客户点
        %         end
        %         if segforlength(i+1) == 0
        %             dot2 = 21;  % 配送中心
        %         else
        %             dot2 = segforlength(i+1);  % 商家或客户点
        %         end 
        %         % dots=[dotss(dot1,:),dotss(dot2,:),];
        %         % if d_matrix(dot1,dot2)==0
        %         % [f_length]=runPlatemo_forOrder_PRMO(@planner_simple_maxhjd_newobj,@simplePRGOCEA,100,30000,20,3,dots,dot1,dot2,1);
        %         % d_matrix(dot1,dot2)=f_length;
        %         % d_matrix(dot2,dot1)=f_length;
        %         % d_hor=f_length;
        %         % else
        %         d_hor=d_matrix(dot1,dot2);
        %         % end
        % 
        %         % 上升（从所有点起飞时都需要上升，因为到达时都会降落）
        %         t_up = h_up_down / v_u;
        %         % 此处写入上升功率计算公式
        %         UAV_w = (UAV_m+load_kg(i))*g; %无人机自重w（T）
        %         v_0 = UAV_w/1.23235;
        %         P_up_i = 79.85628 + UAV_w*(v_u/2+sqrt((v_u/2)^2+v_0));
        %         E_up = P_up_i * t_up;
        % 
        %         % 水平飞行
        %         t_hor = d_hor / v_h;%此处d_hor为航迹规划算法输出的值或矩阵内存调用
        %         %此处写入水平功率计算公式
        %         P_hor_i = 79.85628*(1+0.00020833*v_h^4) + (1.1*UAV_w^(3/2)/sqrt(1.23235))*(sqrt((1+v_h^4)/(4*v_0^4))-v_h^2/(2*v_0))^(1/2) + 0.009249*v_h;
        %         E_hor = P_hor_i * t_hor;
        % 
        %         % 下降（到达所有点时都需要下降，包括depot、商家和客户点）
        %         t_down = h_up_down / v_d;
        %         %此处写入下降功率计算公式
        %         v_d_0=v_d/v_0;
        %         P_down_i = 79.85625 + UAV_w*v_0*(0.974 -1.125*v_d_0 -1.372*v_d_0^2 -1.718*v_d_0^3 -0.655*v_d_0^4);
        %         E_down = P_down_i * t_down;
        % 
        %        % 计算去往每个点路上所花费的时间（单段）
        %        T_k_single = t_up + t_hor + t_down;
        %        % 累积到达时间（从出发到当前“到达点”的总时间）
        %        time_cum = time_cum + T_k_single;
        % 
        %        % ---- 基于"到达点"的时间窗计算时间成本 ----
        %        % 当前段是 segforlength(i) -> segforlength(i+1)
        %        % 到达点索引：0=配送中心(21), 1-10=商家(1-10), 11-20=客户点(11-20)
        %        dest_idx = segforlength(i+1);
        %        if dest_idx == 0
        %            % 配送中心对应 TW_sec 的第21行
        %            tw_row = 21;
        %        elseif dest_idx >= 1 && dest_idx <= n
        %            % 商家（1-10）对应 TW_sec 的第1-10行
        %            tw_row = dest_idx;
        %        elseif dest_idx > n && dest_idx <= 2*n
        %            % 客户点（11-20）对应 TW_sec 的第11-20行
        %            tw_row = dest_idx;  % 客户点11对应第11行，客户点20对应第20行
        %        else
        %            tw_row = 21;  % 默认使用配送中心时间窗
        %        end
        %        aa = TW_sec(tw_row,1);   % 时间窗下界
        %        bb = TW_sec(tw_row,2);   % 时间窗上界
        % 
        %        if time_cum < aa         % 提前到达
        %            T_k = lam1*(aa - time_cum);
        %        elseif time_cum <= bb    % 在时间窗内到达
        %            T_k = time_cum - aa;
        %        else                     % 晚到达
        %            T_k = lam2*(time_cum - bb);
        %        end
        % 
        %        % 累加
        %         % 调试：检查累加前的值
        %         if ~isfinite(T_k) || T_k < 0
        %             debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
        %         end
        %         if ~isfinite(E_up) || ~isfinite(E_hor) || ~isfinite(E_down)
        %             debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
        %         end
        % 
        %         T_k_all = T_k_all + T_k;
        %         E_k = E_k + E_up + E_hor + E_down;
        % 
        %         % 调试：检查累加后的值
        %         if ~isfinite(T_k_all) || ~isfinite(E_k)
        %             debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
        %             debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
        %         end
        %     end
        % 
        %     T_all(k) = T_k_all;%每架飞机的
        %     E_all(k) = E_k;
        % 
        %     % 调试：检查最终的能量和时间值
        %     if ~isfinite(T_all(k)) || T_all(k) < 0
        %         debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
        %     end
        %     if ~isfinite(E_all(k)) || E_all(k) < 0
        %         debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
        %     end
        % end
        % 
        %  % 目标一：总能耗（增加区分度）
        %  baseE_sum = sum(E_all);
        %  if length(E_all) > 1
        %      stdE = real(std(E_all));
        %      rangeE = real(range(E_all));
        %      PopObj(j,1) = baseE_sum + 0.5 * max(0, stdE) + 0.05 * max(0, rangeE) + ...
        %          0.01 * sum(E_all.^2) / length(E_all);
        %  else
        %      PopObj(j,1) = baseE_sum;
        %  end
        %  PopObj(j,1) = real(PopObj(j,1));
        % 
        %  % 目标二：总时间（增加区分度）
        %  baseT_sum = sum(T_all);
        %  if length(T_all) > 1
        %      stdT = real(std(T_all));
        %      rangeT = real(range(T_all));
        %      PopObj(j,2) = baseT_sum + 0.5 * max(0, stdT) + 0.05 * max(0, rangeT) + ...
        %          0.01 * sum(T_all.^2) / length(T_all);
        %  else
        %      PopObj(j,2) = baseT_sum;
        %  end
        %  PopObj(j,2) = real(PopObj(j,2));
        % 
        % % 约束--载重
        % for kk = 1:m
        %     load_k = load_all(kk);
        %     if load_k > maxload
        %         Penalty1=Penalty1+(load_k-maxload);
        %     end
        % end
        % 
        % % 约束--电量
        % for kk = 1:m
        %     power_k = E_all(kk);
        %     if  power_k > maxEC
        %         Penalty2=Penalty2+((power_k-maxEC)/1000);%转为KJ，对应载重约束kg
        %     end
        % 
        % end
        % 
        % % 调试输出：打印前几个解的目标值和惩罚值
        % if j <= 3  % 只打印前3个解的信息
        %     fprintf('解 %d: E_all=[%.2f,%.2f,%.2f], T_all=[%.2f,%.2f,%.2f]\n', ...
        %         j, E_all(1), E_all(2), E_all(3), T_all(1), T_all(2), T_all(3));
        %     fprintf('解 %d: 原始目标值: E=%.2f, T=%.2f\n', j, PopObj(j,1), PopObj(j,2));
        %     fprintf('解 %d: 惩罚值: Penalty1=%.2f, Penalty2=%.2f, Penalty3=%.2f\n', j, Penalty1, Penalty2, Penalty3);
        % end
        % 
        % % 归一化惩罚项，避免目标值过大
        % % 使用相对惩罚，基于目标值的量级
        % baseE = max(1, PopObj(j,1));  % 避免除以0
        % baseT = max(1, PopObj(j,2));  % 避免除以0
        % 
        % % 惩罚系数：根据目标值量级调整（降低惩罚强度，避免过度惩罚）
        % penaltyCoeff = max(1, min(baseE, baseT) / 10000);  % 从1000改为10000，降低惩罚强度
        % 
        % PopObj(j,1)=PopObj(j,1)+penaltyCoeff*(Penalty1+Penalty2)+Penalty3;  % 添加约束惩罚
        % PopObj(j,2)=PopObj(j,2)+penaltyCoeff*(Penalty1+Penalty2)+Penalty3;  % 添加约束惩罚
        % 
        % % 确保目标值为有限值（放宽条件，只在真正无效时才设为1e10）
        % if ~isfinite(PopObj(j,1)) || PopObj(j,1) < 0 || isnan(PopObj(j,1))
        %     PopObj(j,1) = 1e10;  % 设置为一个大的有限值
        % end
        % if ~isfinite(PopObj(j,2)) || PopObj(j,2) < 0 || isnan(PopObj(j,2))
        %     PopObj(j,2) = 1e10;  % 设置为一个大的有限值
        % end
        end
        end


        %% Generate points on the Pareto front
        % function R = GetOptimum(obj,N)
        % % 
        %          R=load('SMEAGLT5.pf');
        % % 
        % end

        %% 诊断方法：返回每架无人机的能耗/用时/惩罚项，用于结果分析
        function Details = CalDetails(obj, PopDec)
        % 返回矩阵 Details（N 行），列为：
        %   1       : Penalty1（载重超限，J当量）
        %   2       : Penalty2（电量超限，J）
        %   3..2+m  : 每架无人机能耗 E_all(1..m)，单位 J
        %   3+m..2+2m: 每架无人机时间代价 T_all(1..m)，单位 s
        %   2+2m+1  : 是否可行（1=可行，0=违约）
        % ---- 参数（与 CalObj 保持一致）----
        N = size(PopDec, 1);
        n = obj.n;  m = obj.m;
        UAV_m=2; demand_q=[0.5;0.7;0.65;0.9;0.35;0.2;0.9;0.1;0.45;0.77;0];
        v_h=10; v_u=5; v_d=3; g=9.8; h_up_down=40;
        rho_air=1.225; A_rotor=0.503;
        Omega=300; R_rotor=0.4; delta_drag=0.012; kappa=1.1; S_FP=0.0151;
        b_blades=4; c_chord=0.0157;
        s_solidity=b_blades*c_chord*R_rotor/(pi*R_rotor^2);
        P_0=(delta_drag/8)*rho_air*s_solidity*A_rotor*Omega^3*R_rotor^3;
        maxload=3; maxEC=800000; loadScale=1e4;  % 与 CalObj 中 maxEC 保持同步
        lam1=0.5; lam2=1.0;
        % 加载地图数据
        current_file_dir=fileparts(mfilename('fullpath'));
        platemo_root=fileparts(fileparts(fileparts(current_file_dir)));
        data_file_path=fullfile(platemo_root,'forOrderNew26','Order_Map','order_data.m');
        if exist(data_file_path,'file'), run(data_file_path); end
        if obj.use_map_penalty
            d_matrix=calculate_distance_matrix_with_obstacles(...
                dotss,forbidden_zones,crowded_zones,noise_zones,...
                penalty_forbidden,penalty_crowded,penalty_noise);
        else
            d_matrix=calculate_distance_matrix_with_obstacles(...
                dotss,[],[],[],0,0,0);
        end
        % 与 CalObj 中 TW_sec 保持同步（原始设置）
        TW_sec=[0,15*60;0,10*60;0,10*60;0,10*60;0,10*60;0,15*60;0,12*60;0,10*60;0,10*60;0,10*60;
                20*60,40*60;10*60,30*60;5*60,25*60;0*60,20*60;5*60,25*60;20*60,40*60;15*60,35*60;10*60,30*60;0*60,20*60;5*60,25*60;
                0,1e11];
        nCols = 2 + 2*m + 1;
        Details = zeros(N, nCols);
        for j = 1:N
            route = PopDec(j,:);
            P1=0; P2=0;
            idx0=find(route==0);
            idx0=[0,idx0,numel(route)+1];
            routes=cell(m,1);
            for k=1:m
                if k<=length(idx0)-1
                    routes{k}=route(idx0(k)+1:idx0(k+1)-1);
                else, routes{k}=[];
                end
            end
            T_all=zeros(1,m); E_all=zeros(1,m); load_all=zeros(1,m);
            for k=1:m
                seg=routes{k};
                if isempty(seg), continue; end
                segforlength=[0,seg,0];
                npts=length(seg)+2;
                load_at_point=zeros(length(segforlength),1);
                for i=2:length(segforlength)-1
                    pid=segforlength(i);
                    if pid>=1&&pid<=n, load_at_point(i)=load_at_point(i-1)+demand_q(pid);
                    elseif pid>n&&pid<=2*n, load_at_point(i)=load_at_point(i-1)-demand_q(pid-n);
                    else, load_at_point(i)=load_at_point(i-1);
                    end
                end
                load_at_point(end)=0;
                load_kg=load_at_point(1:end-1);
                load_all(k)=max(load_at_point);
                T_k=0; E_k=0; time_cum=0;
                for i=1:npts-1
                    dot1=segforlength(i); if dot1==0, dot1=21; end
                    dot2=segforlength(i+1); if dot2==0, dot2=21; end
                    d_hor=d_matrix(dot1,dot2);
                    UAV_mass=UAV_m+load_kg(i); UAV_w=UAV_mass*g;
                    v_0=sqrt((UAV_mass*g)/(2*rho_air*A_rotor));
                    % 上升
                    t_up=h_up_down/v_u;
                    sqrt_term_up=sqrt((v_u/2)^2+v_0^2);
                    P_up=P_0+UAV_w*(v_u/2+sqrt_term_up);
                    E_up=P_up*t_up;
                    % 水平
                    t_hor=d_hor/v_h;
                    mu=v_h/(Omega*R_rotor);
                    P_pro_h=P_0*(1+3*mu^2);
                    sqrt_inner=sqrt(1+v_h^4/(4*v_0^4))-v_h^2/(2*v_0^2);
                    sqrt_term_hor=sqrt(max(0,sqrt_inner));
                    P_i_h=kappa*UAV_w*v_0*sqrt_term_hor;
                    P_par=0.5*rho_air*S_FP*v_h^3;
                    E_hor=(P_pro_h+P_i_h+P_par)*t_hor;
                    % 下降
                    t_down=h_up_down/v_d;
                    v_d_0=v_d/v_0;
                    P_down=P_0+UAV_w*v_0*(0.974-1.125*v_d_0-1.372*v_d_0^2-1.718*v_d_0^3-0.655*v_d_0^4);
                    E_down=max(0,P_down*t_down);
                    T_k_single=t_up+t_hor+t_down;
                    time_cum=time_cum+T_k_single;
                    dest_idx=segforlength(i+1);
                    if dest_idx==0, tw_row=21;
                    elseif dest_idx>=1&&dest_idx<=n, tw_row=dest_idx;
                    elseif dest_idx>n&&dest_idx<=2*n, tw_row=dest_idx;
                    else, tw_row=21;
                    end
                    aa=TW_sec(tw_row,1); bb=TW_sec(tw_row,2);
                    if time_cum<aa, Tk=lam1*(aa-time_cum);
                    elseif time_cum<=bb, Tk=time_cum-aa;
                    else, Tk=lam2*(time_cum-bb);
                    end
                    T_k=T_k+Tk;
                    E_k=E_k+real(E_up)+real(E_hor)+real(E_down);
                end
                T_all(k)=real(T_k); E_all(k)=real(E_k);
            end
            % 约束
            for kk=1:m
                if load_all(kk)>maxload, P1=P1+(load_all(kk)-maxload)*loadScale; end
                if E_all(kk)>maxEC, P2=P2+(E_all(kk)-maxEC); end
            end
            feasible = double(P1==0 && P2==0);
            Details(j,:) = [P1, P2, E_all, T_all, feasible];
        end
        end  % CalDetails

     
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


