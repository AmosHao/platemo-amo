classdef PreDEMO < ALGORITHM
    % <multi/many> <real/binary/permutation>    
    %------------------------------- Reference --------------------------------
    %V. Palakonda and J. -M. Kang, "Pre-DEMO: Preference-Inspired 
    % Differential Evolution for Multi/Many-Objective Optimization," in 
    % IEEE Transactions on Systems, Man, and Cybernetics: Systems, 
    % vol. 53, no. 12, pp. 7618-7630, Dec. 2023
    %--------------------------------------------------------------------------
    
    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population          = Problem.Initialization();
            [FrontNo,~]    = NDSort(Population.objs,Problem.N); 
            %% Optimization
            while Algorithm.NotTerminated(Population)  
               MatingPool1                     = randi(Problem.N,1,Problem.N);
               MatingPool2                     = randi(Problem.N,1,Problem.N);
               Parent2                         = Population(MatingPool1);
               Parent3                         = Population(MatingPool2);               
               Offspring                       = Operator(Population,Parent2,Parent3,FrontNo,Problem.FE,Problem.maxFE,Problem);
               [Population,FrontNo]            = EnvironmentalSelection(Population,Offspring,Problem.N);              
            end
        end
    end
end