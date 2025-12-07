classdef LZ6_ < PROBLEM
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
            x=PopDec;
            n=obj.D;
            M=obj.M;
            PopObj = zeros(N,M);
            i1 = 4:3:n;
            i2 = 5:3:n;
            i3 = 3:3:n;
            ss=zeros(N,n);
        if size(PopDec,1)>1
             for i = 1:N
                 for j=1:n
                 ss(i,j) = x(i,j) - 2.0*x(i,2)*sin(2.0*pi*x(i,1)+j*pi/n);
                 end
               PopObj(i,1)= cos(0.5*pi*x(i,1))*cos(0.5*pi*x(i,2)) + 2.0/length(i1)*sum(ss(i1).^2.0);
               PopObj(i,2)= cos(0.5*pi*x(i,1))*sin(0.5*pi*x(i,2)) + 2.0/length(i2)*sum(ss(i2).^2.0);
               PopObj(i,3)= sin(0.5*pi*x(i,1))                  + 2.0/length(i3)*sum(ss(i3).^2.0);
             end
        else
               for j=1:n
               ss(1,j) = x(1,j) - 2.0*x(1,2)*sin(2.0*pi*x(1,1)+j*pi/n);
               end
               PopObj(1,1)= cos(0.5*pi*x(1,1))*cos(0.5*pi*x(1,2)) + 2.0/length(i1)*sum(ss(i1).^2.0);
               PopObj(1,2)= cos(0.5*pi*x(1,1))*sin(0.5*pi*x(1,2)) + 2.0/length(i2)*sum(ss(i2).^2.0);
               PopObj(1,3)= sin(0.5*pi*x(1,1))                  + 2.0/length(i3)*sum(ss(i3).^2.0);
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
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
     
 end
 end


