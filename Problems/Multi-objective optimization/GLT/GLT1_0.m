classdef GLT1_0 < PROBLEM
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
        N=obj.N;
        D=obj.D;
        g = zeros(N,1);
        PopObj = zeros(N,2);
        if size(PopDec,1)>1
            for i = 1:size(PopDec,1)
            %g(i,1) = sum((PopDec(i,2:D)-cos(2*pi*PopDec(i,1)+(2:D)*pi/D)).^2);
            g(i,1) = sum((PopDec(i,2:D)-sin(2*pi*PopDec(i,1)+(2:D)*pi/D)).^2);
            PopObj(i,1)= (1+g(i,1))*PopDec(i,1);
            PopObj(i,2)= (1+g(i,1))*(2-PopDec(i,1)-sign(cos(2*pi*PopDec(i,1))));
            
             end
        else
            %g(1,1) = sum((PopDec(1,2:D)-cos(2*pi*PopDec(1,1)+(2:D)*pi/D)).^2);
            g(1,1) = sum((PopDec(1,2:D)-sin(2*pi*PopDec(1,1)+(2:D)*pi/D)).^2);
            PopObj(1,1)= (1+g(1,1))*PopDec(1,1);
            PopObj(1,2)= (1+g(1,1))*(2-PopDec(1,1)-sign(cos(2*pi*PopDec(1,1))));

        end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
                 R=load('GLT1.pf');
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

