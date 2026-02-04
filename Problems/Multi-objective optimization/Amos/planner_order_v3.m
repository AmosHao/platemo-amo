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
                % 增加多样性：使用不同的策略生成初始解
                sequence = [];
                strategy = mod(i, 4);  % 使用4种不同策略循环
                
                for p_idx = 1:length(pair_order)
                    pair = pairs{pair_order(p_idx)};
                    
                    % 策略0-1：商家在前（满足约束）
                    % 策略2：随机顺序（增加多样性，后续会修复）
                    % 策略3：客户在前（违反约束，但增加探索空间）
                    if strategy == 0 || strategy == 1
                        % 策略0-1：商家在前（满足约束）
                        sequence = [sequence, pair(1), pair(2)];  % 商家，客户
                    elseif strategy == 2
                        % 策略2：随机顺序
                        if rand < 0.7  % 70%概率商家在前
                            sequence = [sequence, pair(1), pair(2)];
                        else
                            sequence = [sequence, pair(2), pair(1)];
                        end
                    else
                        % 策略3：客户在前（增加探索空间）
                        sequence = [sequence, pair(2), pair(1)];  % 客户，商家（会违反约束，但后续会修复）
                    end
                end
                
                % 随机选择 m-1 个位置插入0分隔符
                % 位置范围：[1, length(sequence)-1]，确保0不在首尾
                if m > 1 && length(sequence) > 1
                    max_pos = min(length(sequence) - 1, 2*n - 1);
                    if max_pos >= m - 1
                        positions = randperm(max_pos, m - 1);  % 随机选择 m-1 个位置
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
                        % 如果位置不足，在末尾补充0
                        result = sequence;
                        result = [result, zeros(1, m - 1)];
                    end
                else
                    result = sequence;
                end
                
                % 确保长度正确（关键：不能补充0，因为0分隔符数量已经正确）
                if length(result) ~= D
                    if length(result) < D
                        % 如果太短，补充非0值（重复已有的点），而不是补充0
                        nonZeroVals = result(result ~= 0);
                        if ~isempty(nonZeroVals)
                            numNeeded = D - length(result);
                            fillValues = nonZeroVals(randi(length(nonZeroVals), 1, numNeeded));
                            result = [result, fillValues];
                        else
                            % 如果没有非0值，使用1-20的随机值填充
                            fillValues = randi([1, 2*n], 1, D - length(result));
                            result = [result, fillValues];
                        end
                    else
                        % 如果太长，截断（但要确保保留m-1个0分隔符）
                        zeroIdx = find(result == 0);
                        if length(zeroIdx) >= m - 1
                            % 保留前m-1个0分隔符
                            lastZeroIdx = zeroIdx(m-1);
                            if lastZeroIdx < D
                                result = result(1:D);
                            else
                                % 如果前m-1个0的位置超过D，需要重新调整
                                result = result(1:D);
                                % 重新确保0分隔符数量
                                zeroCount = sum(result == 0);
                                if zeroCount < m - 1
                                    numNeeded = (m - 1) - zeroCount;
                                    nonZeroIdx = find(result ~= 0);
                                    if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
                                        insertPos = round(linspace(1, length(nonZeroIdx)-1, numNeeded+1));
                                        insertPos = insertPos(2:end);
                                        for i = length(insertPos):-1:1
                                            pos = nonZeroIdx(insertPos(i));
                                            result = [result(1:pos), 0, result(pos+1:end)];
                                        end
                                    end
                                end
                            end
                        else
                            result = result(1:D);
                        end
                    end
                end
                
                % 最终验证：确保0分隔符数量正确
                zeroCount = sum(result == 0);
                if zeroCount ~= m - 1
                    % 强制修正
                    zeroIdx = find(result == 0);
                    if zeroCount < m - 1
                        numNeeded = (m - 1) - zeroCount;
                        nonZeroIdx = find(result ~= 0);
                        if ~isempty(nonZeroIdx) && length(nonZeroIdx) > 1
                            insertPos = round(linspace(1, length(nonZeroIdx)-1, numNeeded+1));
                            insertPos = insertPos(2:end);
                            for i = length(insertPos):-1:1
                                pos = nonZeroIdx(insertPos(i));
                                result = [result(1:pos), 0, result(pos+1:end)];
                            end
                        else
                            result = [result, zeros(1, numNeeded)];
                        end
                        % 如果插入后长度超过D，截断
                        if length(result) > D
                            result = result(1:D);
                        end
                    elseif zeroCount > m - 1
                        result(zeroIdx(m:end)) = [];
                    end
                end
                
                PopDec(i, :) = result;
            end
            
            % 评估初始种群
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
        % 确保有 m-1 个0分隔符
        if length(idx0) < m - 1
            % 如果0的数量不足，在末尾补充0
            numZerosNeeded = (m - 1) - length(idx0);
            % 找到非0的位置，在适当位置插入0
            nonZeroIdx = find(route ~= 0);
            if ~isempty(nonZeroIdx)
                % 在非0元素之间插入0
                for z = 1:numZerosNeeded
                    if length(nonZeroIdx) > 1
                        insertPos = nonZeroIdx(randi(length(nonZeroIdx)-1)) + 1;
                        route = [route(1:insertPos-1), 0, route(insertPos:end)];
                        nonZeroIdx = find(route ~= 0);
                    else
                        route = [route, 0];
                    end
                end
            else
                % 如果全是0，补充到m-1个
                route = [route, zeros(1, numZerosNeeded)];
            end
            idx0 = find(route == 0);
        elseif length(idx0) > m - 1
            % 如果0的数量过多，只保留前m-1个
            idx0 = idx0(1:m-1);
        end
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
                
                % 调试：检查距离矩阵索引
                if dot1 < 1 || dot1 > size(d_matrix,1) || dot2 < 1 || dot2 > size(d_matrix,2)
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                    if length(debug_stats.invalid_examples) < 5
                        debug_stats.invalid_examples{end+1} = struct('type', 'invalid_dot_index', ...
                            'route', route, 'dot1', dot1, 'dot2', dot2, 'matrixSize', size(d_matrix), ...
                            'seg', seg, 'i', i, 'k', k);
                    end
                    d_hor = 0;  % 使用默认值避免崩溃
                end
                
                % 调试：检查距离值
                if ~isfinite(d_hor) || d_hor < 0
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                    d_hor = max(0, d_hor);  % 修正为0
                end

                % 上升（仅在从 depot 或客户点起飞时）
                t_up = h_up_down / v_u;
                
                % 调试：检查时间计算
                if ~isfinite(t_up) || t_up < 0 || v_u <= 0
                    debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
                    t_up = max(0.1, t_up);  % 使用最小值
                end
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
                % 调试：检查累加前的值
                if ~isfinite(T_k) || T_k < 0
                    debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
                end
                if ~isfinite(E_up) || ~isfinite(E_hor) || ~isfinite(E_down)
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                end
                
                T_k_all = T_k_all + T_k;
                E_k = E_k + E_up + E_hor + E_down;
                
                % 调试：检查累加后的值
                if ~isfinite(T_k_all) || ~isfinite(E_k)
                    debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                end
            end
        
            T_all(k) = T_k_all;%每架飞机的
            E_all(k) = E_k;
            
            % 调试：检查最终的能量和时间值
            if ~isfinite(T_all(k)) || T_all(k) < 0
                debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
            end
            if ~isfinite(E_all(k)) || E_all(k) < 0
                debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
            end
        end
        
        % 目标一：总能耗（增加区分度）
        baseE_sum = sum(E_all);
        if length(E_all) > 1
            % 添加区分度项：标准差、范围、平方和
            stdE = real(std(E_all));
            rangeE = real(range(E_all));
            PopObj(j,1) = baseE_sum + 0.5 * max(0, stdE) + 0.05 * max(0, rangeE) + ...
                0.01 * sum(E_all.^2) / length(E_all);
        else
            PopObj(j,1) = baseE_sum;
        end
        PopObj(j,1) = real(PopObj(j,1));
        
        % 目标二：总时间（增加区分度）
        baseT_sum = sum(T_all);
        if length(T_all) > 1
            % 添加区分度项：标准差、范围、平方和
            stdT = real(std(T_all));
            rangeT = real(range(T_all));
            PopObj(j,2) = baseT_sum + 0.5 * max(0, stdT) + 0.05 * max(0, rangeT) + ...
                0.01 * sum(T_all.^2) / length(T_all);
        else
            PopObj(j,2) = baseT_sum;
        end
        PopObj(j,2) = real(PopObj(j,2)); 
         
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
        
        % 归一化惩罚项，避免目标值过大
        % 使用相对惩罚，基于目标值的量级
        baseE = max(1, sum(E_all));  % 避免除以0
        baseT = max(1, sum(T_all));  % 避免除以0
        
        % 惩罚系数：根据目标值量级调整
        penaltyCoeff = max(1, min(baseE, baseT) / 10000);  % 从1000改为10000，降低惩罚强度
        
        PopObj(j,1)=PopObj(j,1)+penaltyCoeff*(Penalty1+Penalty2)+Penalty3;  % 添加约束惩罚
        PopObj(j,2)=PopObj(j,2)+penaltyCoeff*(Penalty1+Penalty2)+Penalty3;  % 添加约束惩罚
        
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
        if j == N
            fprintf('\n========== 调试统计 [批次完成] ==========\n');
            fprintf('总解数: %d\n', N);
            fprintf('无效解总数: %d (%.1f%%)\n', debug_stats.total_invalid, 100*debug_stats.total_invalid/N);
            fprintf('\n无效原因统计:\n');
            fprintf('  - Inf/NaN: %d\n', debug_stats.invalid_reasons.inf_nan);
            fprintf('  - 负数: %d\n', debug_stats.invalid_reasons.negative);
            fprintf('  - 缺少必需的点: %d\n', debug_stats.invalid_reasons.missing_points);
            fprintf('  - 0分隔符数量错误: %d\n', debug_stats.invalid_reasons.wrong_zeros);
            fprintf('  - 空段: %d\n', debug_stats.invalid_reasons.empty_segments);
            fprintf('  - 能量计算错误: %d\n', debug_stats.invalid_reasons.energy_error);
            fprintf('  - 时间计算错误: %d\n', debug_stats.invalid_reasons.time_error);
            
            if ~isempty(debug_stats.invalid_examples)
                fprintf('\n无效解示例（前%d个）:\n', min(3, length(debug_stats.invalid_examples)));
                for ex_idx = 1:min(3, length(debug_stats.invalid_examples))
                    ex = debug_stats.invalid_examples{ex_idx};
                    fprintf('  示例%d - 类型: %s\n', ex_idx, ex.type);
                    fprintf('    路径长度: %d, 路径前15个: [%s]\n', length(ex.route), ...
                        num2str(ex.route(1:min(15, length(ex.route)))));
                    if isfield(ex, 'missing')
                        fprintf('    缺失的点: [%s]\n', num2str(ex.missing));
                    end
                    if isfield(ex, 'zeroCount')
                        if isfield(ex, 'expected')
                            fprintf('    0分隔符数量: %d (期望: %d)\n', ex.zeroCount, ex.expected);
                        else
                            fprintf('    0分隔符数量: %d (期望: %d)\n', ex.zeroCount, m-1);
                        end
                    end
                    if isfield(ex, 'E_all')
                        fprintf('    E_all: [%s]\n', num2str(ex.E_all));
                    end
                    if isfield(ex, 'T_all')
                        fprintf('    T_all: [%s]\n', num2str(ex.T_all));
                    end
                    if isfield(ex, 'PopObj_before')
                        fprintf('    原始目标值: f1=%.2f, f2=%.2f\n', ex.PopObj_before(1), ex.PopObj_before(2));
                    end
                    if isfield(ex, 'Penalty1')
                        fprintf('    惩罚值: P1=%.2f, P2=%.2f, P3=%.2f\n', ex.Penalty1, ex.Penalty2, ex.Penalty3);
                    end
                end
            end
            fprintf('==========================================\n\n');
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
        % 确保有 m-1 个0分隔符
        if length(idx0) < m - 1
            % 如果0的数量不足，在末尾补充0
            numZerosNeeded = (m - 1) - length(idx0);
            % 找到非0的位置，在适当位置插入0
            nonZeroIdx = find(route ~= 0);
            if ~isempty(nonZeroIdx)
                % 在非0元素之间插入0
                for z = 1:numZerosNeeded
                    if length(nonZeroIdx) > 1
                        insertPos = nonZeroIdx(randi(length(nonZeroIdx)-1)) + 1;
                        route = [route(1:insertPos-1), 0, route(insertPos:end)];
                        nonZeroIdx = find(route ~= 0);
                    else
                        route = [route, 0];
                    end
                end
            else
                % 如果全是0，补充到m-1个
                route = [route, zeros(1, numZerosNeeded)];
            end
            idx0 = find(route == 0);
        elseif length(idx0) > m - 1
            % 如果0的数量过多，只保留前m-1个
            idx0 = idx0(1:m-1);
        end
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
                % 调试：检查累加前的值
                if ~isfinite(T_k) || T_k < 0
                    debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
                end
                if ~isfinite(E_up) || ~isfinite(E_hor) || ~isfinite(E_down)
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                end
                
                T_k_all = T_k_all + T_k;
                E_k = E_k + E_up + E_hor + E_down;
                
                % 调试：检查累加后的值
                if ~isfinite(T_k_all) || ~isfinite(E_k)
                    debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
                    debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
                end
            end
        
            T_all(k) = T_k_all;%每架飞机的
            E_all(k) = E_k;
            
            % 调试：检查最终的能量和时间值
            if ~isfinite(T_all(k)) || T_all(k) < 0
                debug_stats.invalid_reasons.time_error = debug_stats.invalid_reasons.time_error + 1;
            end
            if ~isfinite(E_all(k)) || E_all(k) < 0
                debug_stats.invalid_reasons.energy_error = debug_stats.invalid_reasons.energy_error + 1;
            end
        end
        
         % 目标一：总能耗（增加区分度）
         baseE_sum = sum(E_all);
         if length(E_all) > 1
             stdE = real(std(E_all));
             rangeE = real(range(E_all));
             PopObj(j,1) = baseE_sum + 0.5 * max(0, stdE) + 0.05 * max(0, rangeE) + ...
                 0.01 * sum(E_all.^2) / length(E_all);
         else
             PopObj(j,1) = baseE_sum;
         end
         PopObj(j,1) = real(PopObj(j,1));
         
         % 目标二：总时间（增加区分度）
         baseT_sum = sum(T_all);
         if length(T_all) > 1
             stdT = real(std(T_all));
             rangeT = real(range(T_all));
             PopObj(j,2) = baseT_sum + 0.5 * max(0, stdT) + 0.05 * max(0, rangeT) + ...
                 0.01 * sum(T_all.^2) / length(T_all);
         else
             PopObj(j,2) = baseT_sum;
         end
         PopObj(j,2) = real(PopObj(j,2));
         
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
        
        % 归一化惩罚项，避免目标值过大
        % 使用相对惩罚，基于目标值的量级
        baseE = max(1, PopObj(j,1));  % 避免除以0
        baseT = max(1, PopObj(j,2));  % 避免除以0
        
        % 惩罚系数：根据目标值量级调整（降低惩罚强度，避免过度惩罚）
        penaltyCoeff = max(1, min(baseE, baseT) / 10000);  % 从1000改为10000，降低惩罚强度
        
        PopObj(j,1)=PopObj(j,1)+penaltyCoeff*(Penalty1+Penalty2)+Penalty3;  % 添加约束惩罚
        PopObj(j,2)=PopObj(j,2)+penaltyCoeff*(Penalty1+Penalty2)+Penalty3;  % 添加约束惩罚
        
        % 确保目标值为有限值（放宽条件，只在真正无效时才设为1e10）
        if ~isfinite(PopObj(j,1)) || PopObj(j,1) < 0 || isnan(PopObj(j,1))
            PopObj(j,1) = 1e10;  % 设置为一个大的有限值
        end
        if ~isfinite(PopObj(j,2)) || PopObj(j,2) < 0 || isnan(PopObj(j,2))
            PopObj(j,2) = 1e10;  % 设置为一个大的有限值
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


