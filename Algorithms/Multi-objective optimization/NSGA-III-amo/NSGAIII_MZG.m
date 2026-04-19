classdef NSGAIII_MZG < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% NSGA-III-MZG: 嵌入增量聚类配偶选择（MatingSelectionAMO）+ OperatorAMO
%
% 参数说明：
% Kmax --- 20 --- 目标空间聚类最大数量（用于跨簇选配偶；增量维护 centroid）
% pc --- 0.5  --- 交叉概率（交叉/变异二选一）
% er --- 0.25 --- explorationRate（变异算子中 T4 大步长概率基值）

    methods
        function main(Algorithm,Problem)
            %% Parameter
            [Kmax,pc,er] = Algorithm.ParameterSet(20,0.5,0.25);

            %% Generate the reference points and random population
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);

            %% Optimization
            state = struct(); % 聚类状态（跨代复用；增量维护 centroid/clustTag）
            elite = struct('f1',[],'f2',[]);
            while Algorithm.NotTerminated(Population)
                % 父代1：直接使用当前种群（不做锦标赛抽样，避免 cons 全0 时退化为随机采样）
                popSize = length(Population);
                p1 = 1:popSize;
                % 父代2：目标空间增量聚类跨簇配偶选择（对 p1 严格索引配对 p2）
                [p1,p2,state] = MatingSelectionAMO(Population,popSize,'Kmax',Kmax,'p1',p1,'state',state);

                Parents1  = Population(p1(:));
                Parents2  = Population(p2(:));
                Offspring = OperatorAMO(Problem,Parents1,'P2',Parents2,'pc',pc,'er',er); % N 个子代（按索引严格配对）

                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);

                % Decision-space deduplication before selection (if still enough solutions)
                Union = DeduplicateByDec([Population,Offspring]);
                if length(Union) >= Problem.N
                    Population = EnvironmentalSelection(Union,Problem.N,Z,Zmin);
                else
                    Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
                end

                % Endpoint elite writeback (zero extra FE)
                % [Population,elite] = EndpointEliteWriteback(Population,elite);
            end
        end
    end
end

