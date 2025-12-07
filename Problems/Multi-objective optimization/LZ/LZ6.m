classdef LZ6 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+4; end

            obj.lower    = -2*ones(1,obj.D);
            obj.upper    = 2*ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            obj.lower(1) = 0; obj.upper(1) = 1;
            obj.lower(2) = 0; obj.upper(2) = 1;% 设置x1、x2的下界为0
        end
 

        
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        g = zeros(D,1);
        PopObj = zeros(N,3);
        io      = 4:3:D;
        ie      = 5:3:D;
        i3      = 3:3:D;
        % ss      = x - x(1).^(0.5+1.5*((1:1:n)'-2.0)/(n-2.0));
        % 
        % y       = zeros(2,1);
        % y(1)    = x(1) + 2.0/length(io)*sum(ss(io).^2.0);
        % y(2)    = 1.0 - sqrt(x(1)) + 2.0/length(ie)*sum(ss(ie).^2.0);
        % y       = y';
        
        if size(PopDec,1)>1
             for i = 1:N
                 for j=4:3:D
                 
 % ss      = x - 2.0*x(2)*sin(2.0*pi*x(1)+(1:1:n)'*pi/n);
  % ss      = x - 2.0*PopDec(i,2)*sin(2.0*pi*PopDec(i,1)+j*pi/D);
                 g(j,1)=PopDec(i,j)- 2.0*PopDec(i,2)*sin(2.0*pi*PopDec(i,1)+j*pi/D);
                
                 end
                 for k=5:3:D
                 g(k,1)=PopDec(i,k)- 2.0*PopDec(i,2)*sin(2.0*pi*PopDec(i,1)+k*pi/D); 
                
                 end
                 for h=3:3:D
                 g(h,1)=PopDec(i,h)- 2.0*PopDec(i,2)*sin(2.0*pi*PopDec(i,1)+h*pi/D);                
                 end
                 
    


    % y(1)    = cos(0.5*pi*x(1))*cos(0.5*pi*x(2)) + 2.0/length(i1)*sum(ss(i1).^2.0);
    % y(1)    = cos(0.5*pi*PopDec(i,1))*cos(0.5*pi*PopDec(i,2)) + 2.0/length(io)*L1;
    % y(2)    = cos(0.5*pi*x(1))*sin(0.5*pi*x(2)) + 2.0/length(i2)*sum(ss(i2).^2.0);
    % y(3)    = sin(0.5*pi*x(1))                  + 2.0/length(i3)*sum(ss(i3).^2.0); 
    % y(3)    = sin(0.5*pi*PopDec(i,1))                  + 2.0/length(i3)*L3;

                 PopObj(i,1)=cos(0.5*pi*PopDec(i,1))*cos(0.5*pi*PopDec(i,2)) + 2.0/length(io)*sum(g(io,1).^2.0);
                 PopObj(i,2)=cos(0.5*pi*PopDec(i,1))*sin(0.5*pi*PopDec(i,2)) + 2.0/length(ie)*sum(g(ie,1).^2.0);
                 PopObj(i,3)=sin(0.5*pi*PopDec(i,1))  + 2.0/length(i3)*sum(g(i3,1).^2.0);
             end
        else
           
                 for j=4:3:D
                 
 % ss      = x - 2.0*x(2)*sin(2.0*pi*x(1)+(1:1:n)'*pi/n);
  % ss      = x - 2.0*PopDec(i,2)*sin(2.0*pi*PopDec(i,1)+j*pi/D);
                 g(j,1)=PopDec(1,j)- 2.0*PopDec(1,2)*sin(2.0*pi*PopDec(1,1)+j*pi/D);
                 end
                 for k=5:3:D
                 g(k,1)=PopDec(1,k)-2.0*PopDec(1,2)*sin(2.0*pi*PopDec(1,1)+k*pi/D); 
                 end
                 for h=3:3:D
                 g(h,1)=PopDec(1,h)-2.0*PopDec(1,2)*sin(2.0*pi*PopDec(1,1)+h*pi/D);                
               
                 end
                 


                PopObj(1,1)=cos(0.5*pi*PopDec(1,1))*cos(0.5*pi*PopDec(1,2)) + 2.0/length(io)*sum(g(io,1).^2.0);
                 PopObj(1,2)=cos(0.5*pi*PopDec(1,1))*sin(0.5*pi*PopDec(1,2)) + 2.0/length(ie)*sum(g(ie,1).^2.0);
                 PopObj(1,3)=sin(0.5*pi*PopDec(1,1))                  + 2.0/length(i3)*sum(g(i3,1).^2.0);

        end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
                 R=load('LZ09_F6.pf');
        end
         %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                 R = obj.GetOptimum(100);
            elseif obj.M == 3
                %ezmesh('cos(0.5*pi*x1).*cos(0.5*pi*x2)','cos(0.5*pi*x1).*sin(0.5*pi*x2)','sin(0.5*pi*x1)',[0,1],[0,1],10)
                x1 = linspace(0, 1, 10);
                x2 = linspace(0, 1, 10);
                
                f1 = cos(0.5*pi*x1)' * cos(0.5*pi*x2);
                f2 = cos(0.5*pi*x1)' * sin(0.5*pi*x2);
                f3 = repmat(sin(0.5*pi*x1)', 1, size(f1,1));
                
                R = {f1, f2, f3};
                % a = linspace(0,1,10)';
                % R = {a*a'/2,a*(1-a')/2,(1-a)*ones(size(a'))/2};
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

