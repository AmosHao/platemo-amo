classdef huigui < PROBLEM
% <single> <real>
% Bent cigar function

%------------------------------- Reference --------------------------------
% C .T. Yue, K. V. Price, P. N. Suganthan, J. J. Liang, M. Z. Ali, B. Y.
% Qu, N. H. Awad, and P. P Biswas, Problem definitions and evaluation
% criteria for the CEC 2020 special session and competition on single
% objective bound constrained numerical optimization, Zhengzhou University,
% China and Nanyang Technological University, Singapore, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % properties
    %     O;      % Optimal decision vector
    %     Mat;	% Rotation matrix
    % end
    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 1; end
            if isempty(obj.D); obj.D = 3; end
            obj.lower(1:1)   = ones(1,obj.D);
            obj.upper(1:1)   = 3*ones(1,obj.D);
            obj.lower(2:3)   = ones(1,obj.D);
            obj.upper(2:3)   = 3*ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            N=size(PopDec,1);
            M=obj.M;
            PopObj =zeros(N,M);
            for j = 1:N 
            PopObj(j,1)=120302.821-60.716*PopDec(j,1)+98814.742*PopDec(j,2)-46646.082*PopDec(j,3)-698.068*PopDec(j,1)*PopDec(j,2)+125.475*PopDec(j,1)*PopDec(j,3)+13839.033*PopDec(j,2)*PopDec(j,3);
            end
        end
    end
end