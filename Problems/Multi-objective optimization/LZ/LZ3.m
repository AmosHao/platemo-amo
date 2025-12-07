classdef LZ3 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = obj.M+4; end

            obj.lower    = -1*ones(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            obj.lower(1) = 0; % 设置x1、x2的下界为0
        end
 

        
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        g = zeros(D,1);
        PopObj = zeros(N,2);
        io      = 3:2:D;
        ie      = 2:2:D;
        L1=0;L2=0;
        % ss      = x - x(1).^(0.5+1.5*((1:1:n)'-2.0)/(n-2.0));
        % 
        % y       = zeros(2,1);
        % y(1)    = x(1) + 2.0/length(io)*sum(ss(io).^2.0);
        % y(2)    = 1.0 - sqrt(x(1)) + 2.0/length(ie)*sum(ss(ie).^2.0);
        % y       = y';
        if size(PopDec,1)>1
             for i = 1:N
                 for j=3:2:D
                 k=j-1;
                 % ss1     = x - 0.8*x(1)*cos(6.0*pi*x(1)+(1:1:n)'*pi/n);
                 % ss1     = x - 0.8*PopDec(i,1)*cos(6.0*pi*PopDec(i,1)+j*pi/D);
                 % ss2= x - 0.8*x(1)*sin(6.0*pi*x(1)+(1:1:n)'*pi/n);
                 g(j,1)=PopDec(i,j)- 0.8*PopDec(i,1)*cos(6.0*pi*PopDec(i,1)+j*pi/D);
                 end
                 for k=2:2:D
                 g(k,1)=PopDec(i,k)- 0.8*PopDec(i,1)*sin(6.0*pi*PopDec(i,1)+k*pi/D); 
                 end
                 PopObj(i,1)= PopDec(i,1) + 2.0/length(io)*sum(g(io,1).^2.0);
                 PopObj(i,2)=1.0 - sqrt(PopDec(i,1)) + 2.0/length(ie)*sum(g(ie,1).^2.0);
             end
        else
            
            for j=3:2:D
                 k=j-1;
                 g(j,1)=PopDec(1,j)- 0.8*PopDec(1,1)*cos(6.0*pi*PopDec(1,1)+j*pi/D);
            end
            for k=2:2:D
                 g(k,1)=PopDec(1,k)- 0.8*PopDec(1,1)*sin(6.0*pi*PopDec(1,1)+k*pi/D);
               
            end
                 PopObj(1,1)= PopDec(1,1) + 2.0/length(io)*sum(g(io,1).^2.0);
                 PopObj(1,2)=1.0 - sqrt(PopDec(1,1)) + 2.0/length(ie)*sum(g(ie,1).^2.0);

        end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
                 R=load('LZ09_F3.pf');
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

