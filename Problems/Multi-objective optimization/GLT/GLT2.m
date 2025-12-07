classdef GLT2 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = obj.M+4; end
            obj.lower    = -1*ones(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            obj.lower(1) = 0; % 设置x1的下界为0
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        g = zeros(N,1);
        PopObj = zeros(N,2);
        if size(PopDec,1)>1
             for i = 1:size(PopDec,1)
                g(i,1) = sum((PopDec(i,2:D)-abs(cos(2*pi*PopDec(i,1)+(2:D)*pi/D))).^2);
                PopObj(i,1)= (1+g(i,1))*(1-cos(pi*PopDec(i,1)/2));
                PopObj(i,2)= (1+g(i,1))*(10-10*sin(pi*PopDec(i,1)/2));

                % g         = sum((x(2:pd)-abs(cos(2*pi*x(1)+(2:pd)*pi/pd))).^2);
                % y(1)      = (1+g)*(1-cos(pi*x(1)/2));
                 % y(2)      = (1+g)*(10-10*sin(pi*x(1)/2));
             end
        else
            g(1,1) = sum((PopDec(1,2:D)-abs(cos(2*pi*PopDec(1,1)+(2:D)*pi/D))).^2);
            PopObj(1,1)= (1+g(1,1))*(1-cos(pi*PopDec(1,1)/2));
            PopObj(1,2)= (1+g(1,1))*(10-10*sin(pi*PopDec(1,1)/2));
        end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
                 R=load('GLT2.pf');
        end
         %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                 R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,1,10)';
                R = {a*a'/2,a*(1-a')/2,(1-a)*ones(size(a'))/2};
            else
                R = [];
            end
    %         ParetoFront = obj.GetOptimum(100);
    % 
    % % 绘制 Pareto 前沿散点图
    % scatter(ParetoFront(:,1), ParetoFront(:,2), 'filled');
    % xlabel('Objective 1');
    % ylabel('Objective 2');
    % title('Pareto Front');
    % 
    % % 返回空值，因为图像已经绘制出来了
    % R = [];
        end
     
 end
 end

