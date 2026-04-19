function Offspring = OperatorAMO(Problem,Parents,varargin)
% OperatorAMO - 面向 AMO 编码的结构保持交叉/变异算子（可替换 OperatorGA）
%
% 用法（与 PlatEMO 常见 Operator* 风格一致）：
%   Offspring = OperatorAMO(Problem, Parents)
%   Offspring = OperatorAMO(Problem, Parents, 'pc',0.5, 'er',0.25)
%   Offspring = OperatorAMO(Problem, Parents, 'P2', ExtraParents)
%
% 输入：
%   Problem : PlatEMO 的 PROBLEM 对象（要求提供 Problem.n, Problem.m, Problem.D）
%   Parents : INDIVIDUAL 数组（长度通常为 N），至少包含 .decs
%   ExtraParents : INDIVIDUAL 数组或决策矩阵（行数需与 Parents 一致）。若提供，则交叉时严格按索引 i 使用 ExtraParents(i,:) 作为第二父代
%
% 参数（可选，键值对）：
%   pc : 交叉概率（默认 0.5），每对父代以 pc 进行交叉，否则对 parent1 变异
%   er : explorationRate（默认 0.25），传给 OrderPreservingMutation 的大步长概率基值
%
% 输出：
%   Offspring : INDIVIDUAL 数组，通过 Problem.Evaluation(OffDec) 生成
%
% 说明：
% - 本算子只负责“交叉/变异 + 段结构修复”，不做目标空间聚类选择（聚类交配请用 MatingSelectionAMO）
% - 需要依赖 A-amos 目录中的：
%   SegmentBasedCrossover.m, OrderPreservingMutation.m, OrderSegmentRepair.m

    % ---- parse args ----
    pc = 0.5;
    er = 0.25;
    P2 = [];
    if ~isempty(varargin)
        for i = 1:2:length(varargin)
            key = lower(string(varargin{i}));
            val = varargin{i+1};
            switch key
                case "pc"
                    pc = val;
                case {"er","explorationrate"}
                    er = val;
                case {"p2","extraparents","parent2","parents2"}
                    P2 = val;
            end
        end
    end

    n      = Problem.n;
    m      = Problem.m;
    varDim = Problem.D;

    ParentDec = Parents.decs;
    if isempty(ParentDec)
        Offspring = Problem.Evaluation(ParentDec);
        return;
    end

    % Optional paired second parents (strict index match)
    if ~isempty(P2)
        if isa(P2(1),'SOLUTION') || isa(P2(1),'INDIVIDUAL')
            P2Dec = P2.decs;
        else
            P2Dec = P2;
        end
        if size(P2Dec,1) ~= size(ParentDec,1)
            error('OperatorAMO:InvalidP2Size','ExtraParents (P2) must have the same number of rows as Parents.');
        end
        if size(P2Dec,2) ~= size(ParentDec,2)
            error('OperatorAMO:InvalidP2Dim','ExtraParents (P2) must have the same decision dimension as Parents.');
        end
    else
        P2Dec = [];
    end

    % 每个父代生成 1 个子代（交叉：使用 P2(i,:) 或在 Parents 内随机挑一个 != i）
    N = size(ParentDec,1);
    if N < 1
        Offspring = Problem.Evaluation(ParentDec);
        return;
    end

    OffDec = zeros(N,varDim);
    for i = 1:N
        p1 = ParentDec(i,:);
        useCrossover = (rand < pc) && (N >= 2);
        if useCrossover
            if ~isempty(P2Dec)
                p2 = P2Dec(i,:);
            else
                idx2 = randi(N-1);
                if idx2 >= i
                    idx2 = idx2 + 1;
                end
                p2 = ParentDec(idx2,:);
            end
            child = SegmentBasedCrossover(p1,p2,n,m,varDim);
        else
            child = OrderPreservingMutation(p1,m,varDim,n,er);
        end
        child = OrderSegmentRepair(child,n,m,varDim);
        if length(child) ~= varDim
            if length(child) < varDim
                child = [child,zeros(1,varDim-length(child))];
            else
                child = child(1:varDim);
            end
        end
        OffDec(i,:) = child;
    end

    Offspring = Problem.Evaluation(OffDec);
end

