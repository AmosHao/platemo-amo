classdef MOEADUR_G < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% MOEA/D with update when required (OperatorAMO variant)
% start  ---  0.2 --- Start adaptation
% finish --- 0.93 --- Finish adaptation
% K      ---   10 --- Number of divisions of the objective space
% pc     ---  0.5 --- Crossover probability in OperatorAMO
% er     --- 0.25 --- Exploration rate in OperatorAMO
%
%------------------------------- Reference --------------------------------
% L. R. de Farias, A. F. Araujo, A decomposition-based many-objective
% evolutionary algorithm updating weights when required, Swarm and
% Evolutionary Computation, 2022, 68: 100980.
%--------------------------------------------------------------------------
	
	methods
		function main(Algorithm,Problem)
            %% Parameter setting
            [start, finish, K, pc, er] = Algorithm.ParameterSet(0.2, 0.93, 10, 0.5, 0.25);

            %% MOEA/D parameters
            delta = 0.9;            % The probability of choosing parents locally
            nr = 2;                 % Maximum number of solutions replaced by each offspring
            T = ceil(Problem.N/10); % Size of neighborhood
            mini_generation = 1;    % The number of generations carried out within the objective space division method

            %% Generate the weight vectors
            [W,Problem.N] = UniformlyRandomlyPoint(Problem.N,Problem.M);
            W = 1./W./repmat(sum(1./W,2),1,size(W,2)); % WS-Transformation on W
            W_URP = W;

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);
            EP         = Population(NDSort(Population.objs,1)==1); % external population

            WhenDoesItStart = floor(start*(Problem.maxFE/Problem.N));
            whenItEnds      = floor(finish*(Problem.maxFE/Problem.N));

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % For each solution
                Offsprings(1:Problem.N) = SOLUTION();
                chosenNeighborhood = rand(Problem.N,1) < delta;
                for i = 1 : Problem.N
                    % Choose the parents (index list P)
                    if chosenNeighborhood(i)==1
                        P = B(i,randperm(size(B,2)));
                    else
                        P = randperm(Problem.N);
                    end

                    % Generate ONE offspring (match OperatorGAhalf behavior)
                    Parent1 = Population(P(1));
                    Parent2 = Population(P(2));
                    child   = OperatorAMO(Problem,Parent1,'P2',Parent2,'pc',pc,'er',er);
                    Offsprings(i) = child;

                    % Update the ideal point
                    Z = min(Z,Offsprings(i).obj);

                    % Update the solutions in P by Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offsprings(i).obj-Z),length(P),1).*W(P,:),[],2);
                    Population(P(find(g_old>=g_new,nr))) = Offsprings(i);
                end

                if Problem.FE/Problem.N == WhenDoesItStart
                    X = unique(Population.objs,'rows');
                    X = X(NDSort(X,1)==1,:);
                    X = normalize(X,'norm');
                    spreading_index = norm(X)/4;

                    fun_threshold=[-1.989e-05;0.0002034;0.03376;0.2373];
                    threshold = polyval(fun_threshold,Problem.M);

                    if spreading_index<=threshold
                        period = 12;
                        nus    = 0.25;
                    else
                        period = 28;
                        nus    = 0.075;
                    end

                    rate_update_weight = round(nus*Problem.N);

                    fun_rho = [-0.4707,0.8644,-0.1508,0.05745];
                    rho     = polyval(fun_rho,spreading_index);

                    I_old = max(abs((Population.objs-repmat(Z,Problem.N,1)).*W),[],2);
                end

                % EP update
                if Problem.FE/Problem.N <= whenItEnds
                    EP = [EP,Offsprings];
                    % Decision-space deduplication for archive (does not change FE)
                    EP = DeduplicateByDec(EP);
                end

                if  Problem.FE/Problem.N >= WhenDoesItStart && Problem.FE/Problem.N <= whenItEnds
                    if ~mod(Problem.FE/Problem.N,period)
                        % Improvement Metric
                        I_new = max(abs((Population.objs-repmat(Z,Problem.N,1)).*W),[],2);
                        improvement_Metric = mean(1-(I_new./I_old));
                        I_old = I_new;

                        if abs(improvement_Metric)<=rho
                            % EP update
                            EP = DeduplicateByDec(EP);
                            EP = unique(EP);
                            EP = EP(NDSort(EP.objs,1)==1);

                            % adaptive weight adjustment
                            [Population,W,B] = updateWeight(Population,W,Z,T,EP,rate_update_weight);
                            I_old = max(abs((Population.objs-repmat(Z,Problem.N,1)).*W),[],2);

                            % division of objective space
                            [Population,Z] = space_divide(Problem,EP,Population,Z,W_URP,W,K,mini_generation);
                            EP = [];
                        end
                    end
                end
            end
        end
    end
end

