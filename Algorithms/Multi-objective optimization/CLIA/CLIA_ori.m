classdef CLIA_ori < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% CLIA baseline with RepairOrderEncoding (dec -> repair -> evaluation)
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
            [stable_threshold, delta] = Algorithm.ParameterSet([0, 0, 0], 2 * Problem.N);
            [Z, P, A, S, SVM] = initialize(Problem);
            while Algorithm.NotTerminated(S)
                MatingPool = TournamentSelection(2, Problem.N, sum(max(0, P.cons), 2));

                ParentsDec = P(MatingPool).decs;
                OffDec     = OperatorGA(Problem,ParentsDec);
                OffDec     = RepairOrderEncoding(Problem,OffDec);
                Offspring  = Problem.Evaluation(OffDec);

                A = update_archive(A, [P, Offspring], Z, ceil(0.33 * Problem.M * Problem.N), Problem);
                % SELECTION OF INDIVIDUALS
                [P, ICA, ICN] = cascade_cluster([P, Offspring], Z, 'PDM', Problem.N, Problem.FE < Problem.maxFE);
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

