classdef CLIA_G < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% CLIA (OperatorAMO variant)
%
% 参数说明：
% pc --- 0.5  --- 交叉概率（交叉/变异二选一）
% er --- 0.25 --- explorationRate（变异算子中 T4 大步长概率基值）
%
%------------------------------- Reference --------------------------------
% H. Ge, M. Zhao, L. Sun, Z. Wang, G. Tan, Q. Zhang, and C. L. P. Chen, A
% many-objective evolutionary algorithm with two interacting processes:
% Cascade clustering and reference point incremental learning, IEEE
% Transactions on Evolutionary Computation, 2019, 23(4): 572-586.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            global stable_threshold delta crowding_pick_flag;
            [stable_threshold, delta, pc, er] = Algorithm.ParameterSet([0, 0, 0], 2 * Problem.N, 0.5, 0.25);
            [Z, P, A, S, SVM] = initialize(Problem);
            while Algorithm.NotTerminated(S)
                MatingPool = TournamentSelection(2, Problem.N, sum(max(0, P.cons), 2));
                Parents    = P(MatingPool);
                Offspring  = OperatorAMO(Problem, Parents, 'pc', pc, 'er', er);

                % Decision-space deduplication before downstream selection/archive
                Union = DeduplicateByDec([P, Offspring]);

                A = update_archive(A, Union, Z, ceil(0.33 * Problem.M * Problem.N), Problem);
                % SELECTION OF INDIVIDUALS
                [P, ICA, ICN] = cascade_cluster(Union, Z, 'PDM', Problem.N, Problem.FE < Problem.maxFE);
                % ADAPTATION OF REFERENCE VECTORS
                [Z, SVM] = incremental_learn(Z, ICA, ICN, A, SVM, Problem);
                if Problem.FE >= Problem.maxFE && crowding_pick_flag
                    S = crowding_pick(update_archive(A, P, Z, [], Problem), Problem.N, 'precise');
                else
                    S = P;
                end
            end
        end
    end
end

