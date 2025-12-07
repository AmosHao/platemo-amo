classdef FLYPRO2 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D =30; end
            D=obj.D;
            n=D/3; 
            obj.lower(1) = 0;
            obj.lower(2:n-1) = zeros(1,n-2);
            obj.lower(n) = 19.9;

            obj.lower(n+1) =0;
            obj.lower(n+2:2*n-1) =zeros(1,n-2);
            obj.lower(2*n) =19.9;

            obj.lower(2*n+1) =0;
            obj.lower(2*n+2:3*n-1) =0.15*ones(1,n-2);
            obj.lower(3*n) =0.299;

            obj.upper(1) = 0.1;
            obj.upper(2:n-1) = 20*ones(1,n-2);
            obj.upper(n) = 20.0;

            obj.upper(n+1) =0.1;
            obj.upper(n+2:2*n-1) =20*ones(1,n-2);
            obj.upper(2*n) =20.0;

            obj.upper(2*n+1) =0.1;
            obj.upper(2*n+2:3*n-1) =0.3*ones(1,n-2);
            obj.upper(3*n) =0.300;

            % obj.lower(1:n) = zeros(1,n);
            % obj.lower(n+1:2*n) =zeros(1,n);
            % obj.lower(2*n+1:3*n) =150*ones(1,n);
            % obj.upper(1:n) = 200*ones(1,n);
            % obj.upper(n+1:2*n) = 200*ones(1,n);
            % obj.upper(2*n+1:3*n) = 300*ones(1,n);
            obj.encoding = ones(1,obj.D);
         
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        M=obj.M;
        n=D/3;
        PopObj =zeros(N,M);
        p1=0;p2=0;p3=0;
        q1=20;q2=20;q3=0.3;

         if size(PopDec,1)>1
             for i = 1:N 
               
PopObj(i,1)=sum(sqrt((PopDec(i,2:n)-PopDec(i,1:n-1)).^2+(PopDec(i,n+2:2*n)-PopDec(i,n+1:2*n-1)).^2+(PopDec(i,2*n+2:3*n)-PopDec(i,2*n+1:3*n-1)).^2+(PopDec(i,1)-p1).^2+(PopDec(i,n+1)-p2).^2+(PopDec(i,2*n+1)-p3).^2+(PopDec(i,n)-q1).^2+(PopDec(i,2*n)-q2).^2+(PopDec(i,3*n)-q3).^2));
g1=1-(((PopDec(i,1)-p1)*(PopDec(i,2)-PopDec(i,1)))+((PopDec(i,n+1)-p2)*(PopDec(i,n+2)-PopDec(i,n+1))))/sqrt(((PopDec(i,1)-p1).^2+(PopDec(i,n+1)-p2).^2)*((PopDec(i,2)-PopDec(i,1)).^2+(PopDec(i,n+2)-PopDec(i,n+1)).^2));
g2=1-(((PopDec(i,n)-PopDec(i,n-1))*(q1-PopDec(i,n)))+((PopDec(i,2*n)-PopDec(i,2*n-1))*(q2-PopDec(i,2*n))))/sqrt(((PopDec(i,n)-PopDec(i,n-1)).^2+(PopDec(i,2*n)-PopDec(i,2*n-1)).^2)*((q1-PopDec(i,n)).^2+(q2-PopDec(i,2*n)).^2));
PopObj(i,2)=sum(PopDec(i,2*n+1:3*n));

PopObj(i,3)=g1+g2+sum(1-(((PopDec(i,2:n-1)-PopDec(i,1:n-2)).*(PopDec(i,3:n)-PopDec(i,2:n-1)))+((PopDec(i,n+2:2*n-1)-PopDec(i,n+1:2*n-2)).*(PopDec(i,n+3:2*n)-PopDec(i,n+2:2*n-1))))/sqrt(((PopDec(i,2:n-1)-PopDec(i,1:n-2)).^2+(PopDec(i,n+2:2*n-1)-PopDec(i,n+1:2*n-2)).^2).*((PopDec(i,3:n)-PopDec(i,2:n-1)).^2+(PopDec(i,n+3:2*n)-PopDec(i,n+2:2*n-1)).^2)));
   
             end
         else
PopObj(1,1)=sum(sqrt((PopDec(1,2:n)-PopDec(1,1:n-1)).^2+(PopDec(1,n+2:2*n)-PopDec(1,n+1:2*n-1)).^2+(PopDec(1,2*n+2:3*n)-PopDec(1,2*n+1:3*n-1)).^2+(PopDec(1,1)-p1).^2+(PopDec(1,n+1)-p2).^2+(PopDec(1,2*n+1)-p3).^2+(PopDec(1,n)-q1).^2+(PopDec(1,2*n)-q2).^2+(PopDec(1,3*n)-q3).^2));
g1=1-(((PopDec(1,1)-p1)*(PopDec(1,2)-PopDec(1,1)))+((PopDec(1,n+1)-p2)*(PopDec(1,n+2)-PopDec(1,n+1))))/sqrt(((PopDec(1,1)-p1).^2+(PopDec(1,n+1)-p2).^2)*((PopDec(1,2)-PopDec(1,1)).^2+(PopDec(1,n+2)-PopDec(1,n+1)).^2));
g2=1-(((PopDec(1,n)-PopDec(1,n-1))*(q1-PopDec(1,n)))+((PopDec(1,2*n)-PopDec(1,2*n-1))*(q2-PopDec(1,2*n))))/sqrt(((PopDec(1,n)-PopDec(1,n-1)).^2+(PopDec(1,2*n)-PopDec(1,2*n-1)).^2)*((q1-PopDec(1,n)).^2+(q2-PopDec(1,2*n)).^2));
PopObj(1,2)=sum(PopDec(1,2*n+1:3*n));

PopObj(1,3)=g1+g2+sum(1-(((PopDec(1,2:n-1)-PopDec(1,1:n-2)).*(PopDec(1,3:n)-PopDec(1,2:n-1)))+((PopDec(1,n+2:2*n-1)-PopDec(1,n+1:2*n-2)).*(PopDec(1,n+3:2*n)-PopDec(1,n+2:2*n-1))))/sqrt(((PopDec(1,2:n-1)-PopDec(1,1:n-2)).^2+(PopDec(1,n+2:2*n-1)-PopDec(1,n+1:2*n-2)).^2).*((PopDec(1,3:n)-PopDec(1,2:n-1)).^2+(PopDec(1,n+3:2*n)-PopDec(1,n+2:2*n-1)).^2)));
        %     g(1,1) = sum((PopDec(1,2:D)-abs(cos(2*pi*PopDec(1,1)+(2:D)*pi/D))).^2);
        %     PopObj(1,1)= (1+g(1,1))*(1-cos(pi*PopDec(1,1)/2));
        %     PopObj(1,2)= (1+g(1,1))*(10-10*sin(pi*PopDec(1,1)/2));
         end
 end     
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        n=D/3;
        PopCon = zeros(N,2);
        p1=0;p2=0;p3=0;
        q1=20;q2=20;q3=0.3;
            % t      = X(:,2) - sin(6*pi*X(:,1)+2*pi/size(X,2)) - 0.5*X(:,1) + 0.25;
            % PopCon = -t./(1+exp(4*abs(t)));
        if size(PopDec,1)>1
        for i = 1:N
            PopCon(i,1)=cos(pi/3)-(((PopDec(i,1)-p1)*(PopDec(i,2)-PopDec(i,1)))+((PopDec(i,n+1)-p2)*(PopDec(i,n+2)-PopDec(i,n+1))))/sqrt(((PopDec(i,1)-p1).^2+(PopDec(i,n+1)-p2).^2)*((PopDec(i,2)-PopDec(i,1)).^2+(PopDec(i,n+2)-PopDec(i,n+1)).^2));
            PopCon(i,1)=cos(pi/3)-(((PopDec(i,n)-PopDec(i,n-1))*(q1-PopDec(i,n)))+((PopDec(i,2*n)-PopDec(i,2*n-1))*(q2-PopDec(i,2*n))))/sqrt(((PopDec(i,n)-PopDec(i,n+1)).^2+(PopDec(i,2*n)-PopDec(i,2*n-1)).^2)*((q1-PopDec(i,n)).^2+(q2-PopDec(i,2*n)).^2));
            PopCon(i,2)=(PopDec(i,2*n+1)-p3)/sqrt((PopDec(i,1)-p1).^2+(PopDec(i,n+1)-p2).^2)-tan(pi/4);
            PopCon(i,2)=(PopDec(i,3*n)-PopDec(i,3*n-1))/sqrt((PopDec(i,n)-PopDec(i,n-1)).^2+(PopDec(i,2*n)-PopDec(i,2*n-1)).^2)-tan(pi/4);
            for j=2:n-1
            PopCon(i,1)=cos(pi/3)-(((PopDec(i,j)-PopDec(i,j-1))*(PopDec(i,j+1)-PopDec(i,j)))+((PopDec(i,n+j)-PopDec(i,n+j-1))*(PopDec(i,n+j+1)-PopDec(i,n+j))))/sqrt(((PopDec(i,j)-PopDec(i,j-1)).^2+(PopDec(i,n+j)-PopDec(i,n+j-1)).^2)*((PopDec(i,j+1)-PopDec(i,j)).^2+(PopDec(i,n+j+1)-PopDec(i,n+j)).^2));
            PopCon(i,2)=(PopDec(i,2*n+j)-PopDec(i,2*n+j-1))/sqrt((PopDec(i,j)-PopDec(i,j-1)).^2+(PopDec(i,n+j)-PopDec(i,n+j-1)).^2)-tan(pi/4);
            end
        end
        else
            PopCon(1,1)=cos(pi/3)-(((PopDec(1,1)-p1)*(PopDec(1,2)-PopDec(1,1)))+((PopDec(1,n+1)-p2)*(PopDec(1,n+2)-PopDec(1,n+1))))/sqrt(((PopDec(1,1)-p1).^2+(PopDec(1,n+1)-p2).^2)*((PopDec(1,2)-PopDec(1,1)).^2+(PopDec(1,n+2)-PopDec(1,n+1)).^2));
            PopCon(1,1)=cos(pi/3)-(((PopDec(1,n)-PopDec(1,n-1))*(q1-PopDec(1,n)))+((PopDec(1,2*n)-PopDec(1,2*n-1))*(q2-PopDec(1,2*n))))/sqrt(((PopDec(1,n)-PopDec(1,n+1)).^2+(PopDec(1,2*n)-PopDec(1,2*n-1)).^2)*((q1-PopDec(1,n)).^2+(q2-PopDec(1,2*n)).^2));
            PopCon(1,2)=(PopDec(1,2*n+1)-p3)/sqrt((PopDec(1,1)-p1).^2+(PopDec(1,n+1)-p2).^2)-tan(pi/4);
            PopCon(1,2)=(PopDec(1,3*n)-PopDec(1,3*n-1))/sqrt((PopDec(1,n)-PopDec(1,n-1)).^2+(PopDec(1,2*n)-PopDec(1,2*n-1)).^2)-tan(pi/4);
            for j=2:n-1
            PopCon(1,1)=cos(pi/3)-(((PopDec(1,j)-PopDec(1,j-1))*(PopDec(1,j+1)-PopDec(1,j)))+((PopDec(1,n+j)-PopDec(1,n+j-1))*(PopDec(1,n+j+1)-PopDec(1,n+j))))/sqrt(((PopDec(1,j)-PopDec(1,j-1)).^2+(PopDec(1,n+j)-PopDec(1,n+j-1)).^2)*((PopDec(1,j+1)-PopDec(1,j)).^2+(PopDec(1,n+j+1)-PopDec(1,n+j)).^2));
            PopCon(1,2)=(PopDec(1,2*n+j)-PopDec(1,2*n+j-1))/sqrt((PopDec(1,j)-PopDec(1,j-1)).^2+(PopDec(1,n+j)-PopDec(1,n+j-1)).^2)-tan(pi/4);
            end
        end     
        end

       %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
                 R=[20,1,0.0001];
                 %R=[5,5,0.01];
        end
    % %% Generate points on the Pareto front
    %     function R = GetOptimum(obj,N)
    %         R = UniformPoint(N,obj.M)/2;
    %     end
    %     %% Generate the image of Pareto front
    %     function R = GetPF(obj)
    %         if obj.M == 2
    %             R = obj.GetOptimum(100);
    %         elseif obj.M == 3
    %             a = linspace(0,1,10)';
    %             R = {a*a'/2,a*(1-a')/2,(1-a)*ones(size(a'))/2};
    %         else
    %             R = [];
    %         end
    %     end
     
 end
end


