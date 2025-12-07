classdef LZ8 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = obj.M+4; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);

        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            N=size(PopDec,1);
            x=PopDec;
            n=obj.D;
            M=obj.M;
            PopObj = zeros(N,M);
            io      = 3:2:n;
            ie      = 2:2:n;
            yy=zeros(N,n);ss1=zeros(N,n);ss2=zeros(N,n);
            if size(PopDec,1)>1
                for i=1:N
                    for j=1:n
                    yy(i,j)= x(i,j) - x(i,1).^(0.5+1.5*(j-2)/(n-2.0));
                    ss1(i,j)= yy(i,j)*yy(i,j);
                    ss2(i,j)= cos(20*pi*yy(i,j)./sqrt(j));
                    end
                PopObj(i,1) = x(i,1) + 2.0/length(io)*(4*sum(ss1(i,io))-2*prod(ss2(i,io))+2);
                PopObj(i,2)  = 1.0 - sqrt(x(i,1))  + 2.0/length(ie)*(4*sum(ss1(i,ie))-2*prod(ss2(i,ie))+2);
                end
            else
                for j=1:n
                yy(1,j)= x(1,j) - x(1,1).^(0.5+1.5*(j-2)/(n-2.0));
                ss1(1,j)= yy(1,j)*yy(1,j);
                ss2(1,j)= cos(20*pi*yy(1,j)./sqrt(j));
                end
                PopObj(1,1) = x(1,1) + 2.0/length(io)*(4*sum(ss1(1,io))-2*prod(ss2(1,io))+2);
                PopObj(1,2)  = 1.0 - sqrt(x(1,1))  + 2.0/length(ie)*(4*sum(ss1(1,ie))-2*prod(ss2(1,ie))+2);
            end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
                 R=load('LZ09_F8.pf');
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

