classdef planner_curve < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D =63; end
            D=obj.D;

            n=(D-3)/6;
            % obj.lower(1:2*n) = zeros(1,20);
            obj.lower(2*n+1:3*n) = 0.01*ones(1,10);
            obj.lower(3*n+1:5*n+2) = zeros(1,22);
            obj.lower(5*n+3:6*n+3) = 0.01*ones(1,11);

            % obj.upper(1:2*n) = 2*ones(1,20);
            obj.upper(2*n+1:3*n) = 0.3*ones(1,10);
            obj.upper(3*n+1:5*n+2) = 2*ones(1,22);
            obj.upper(5*n+3:6*n+3) = 0.3*ones(1,11);
            for cc=1:10
            obj.lower(cc) = (cc-1)*0.2*ones(1);
            obj.lower(cc+10)=(cc-1)*0.2*ones(1);
            obj.lower(cc+20)=(cc-1)*0.03*ones(1);
            end
            for dd=1:10
            obj.upper(dd) = dd*0.2*ones(1);
            obj.upper(dd+10) = dd*0.2*ones(1);
            obj.upper(dd+20) = dd*0.03*ones(1);
            end            


            obj.encoding = ones(1,obj.D);
         
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;

        n=(D-3)/6;
        M=obj.M;
        PopObj =zeros(N,M);

        % 半径惩罚参数
        penaltyLarge = 1000;  % 超过上限的惩罚
        penaltySmall = 1000;  % 低于下限的惩罚
        radiusMin = 0.01;%半径上下限
        radiusMax = 1;

        if size(PopDec,1)>1
        for j = 1:N 
        Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;Penalty7=0;
        xStart=[0,0,0]; xEnd=[2,2,0.3];                              % 起点和终点
        z=PopDec(j,2*n+1:3*n); y=PopDec(j,n+1:2*n); x=PopDec(j,1:n);
        x=[xStart(1),x,xEnd(1)]; y=[xStart(2),y,xEnd(2)]; z=[xStart(3),z,xEnd(3)]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        nNodes=n+2;
        ox=PopDec(j,3*n+1:4*n+1); oy=PopDec(j,4*n+2:5*n+2); oz=PopDec(j,5*n+3:6*n+3);
       
        dis=zeros(1,nNodes-1);
        totalCircumference = 0;
        radii = zeros(1, nNodes-1);

        for i=1:nNodes-1
            if x(i+1)-x(i)<0
                Penalty1=Penalty1+10000*abs((x(i)-x(i+1)));%惩罚项1，保证前进
            end
            if y(i+1)-y(i)<0
                Penalty1=Penalty1+10000*abs((y(i+1)-y(i)));
            end
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)<0.1
                Penalty2=Penalty2+dis(i);%惩罚项2，限制两点之间距离大于0.1(2)
            end

            % 计算相邻航迹点到圆心的距离
            dist1 = sqrt((x(i) - ox(i))^2 + (y(i) - oy(i))^2 + (z(i) - oz(i))^2);
            dist2 = sqrt((x(i+1) - ox(i))^2 + (y(i+1) - oy(i))^2 + (z(i+1) - oz(i))^2);
            
            % 计算相邻航迹点之间的距离
            dis3 = sqrt((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2 + (z(i+1) - z(i))^2);
            
            % 惩罚项3：确保相邻两点到圆心的距离相等
            Penalty3 = Penalty3 + 10000*abs(dist1 - dist2);
            
            % 惩罚项4：确保两点之间的距离小于半径的二倍
            radius = dist1; % 半径
            if dis3 >= 2 * radius
                Penalty4 = Penalty4 + 10*(dis3 - 2 * radius); % 超过部分的惩罚
            end

            % 计算圆的周长
            circumference = 2 * pi * radius;
            
            % 累加周长
            totalCircumference = totalCircumference + circumference;

            radii(i) = radius;

            % 惩罚项5
            if radius < radiusMin
                Penalty5 = Penalty5 + penaltySmall * (radiusMin - radius);
            elseif radius > radiusMax
                Penalty5 = Penalty5 + penaltyLarge * (radius - radiusMax);
            end
        end

        PopObj(j,1)=totalCircumference; %目标值1，圆的周长之和   
        
        % 计算半径均值
        meanRadius = mean(radii);
        % 计算均值的倒数
        if meanRadius > 0
            reciprocalMeanRadius = 1 / meanRadius;
        else
            reciprocalMeanRadius = Inf; % 如果均值为0，倒数设为无穷大
        end

        PopObj(j,2) = reciprocalMeanRadius;%目标值1，圆的半径均值倒数 


        if PopObj(j,1)>40
            Penalty6=Penalty6+(PopObj(j,1)-40);%惩罚项6，限制总航迹长度不能大于40（500）
        end

        for i=2:nNodes                                                     % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            if delta_z/norm(a_i)>1
            Penalty7=Penalty7+delta_z/norm(a_i);%惩罚项7，限制爬升/俯冲角tan值小于1
            end
        end

        PopObj(j,1)=PopObj(j,1)+0.07*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6+Penalty7);
        PopObj(j,2)=PopObj(j,2)+0.07*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6+Penalty7);
        
      end
   
            
   else
        Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;Penalty7=0;
        xStart=[0,0,0]; xEnd=[2,2,0.3];                              % 起点和终点
        z=PopDec(1,2*n+1:3*n); y=PopDec(1,n+1:2*n); x=PopDec(1,1:n);
        x=[xStart(1),x,xEnd(1)]; y=[xStart(2),y,xEnd(2)]; z=[xStart(3),z,xEnd(3)]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        nNodes=n+2;
        ox=PopDec(1,3*n+1:4*n+1); oy=PopDec(1,4*n+2:5*n+2); oz=PopDec(1,5*n+3:6*n+3);
       
        dis=zeros(1,nNodes-1);
        totalCircumference = 0;
        radii = zeros(1, nNodes-1);

        for i=1:nNodes-1
            if x(i+1)-x(i)<0
                Penalty1=Penalty1+(x(i)-x(i+1));%惩罚项1，保证前进
            end
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)<0.1
                Penalty2=Penalty2+dis(i);%惩罚项2，限制两点之间距离大于0.1(2)
            end

                % 计算相邻航迹点到圆心的距离
                dist1 = sqrt((x(i) - ox(i))^2 + (y(i) - oy(i))^2 + (z(i) - oz(i))^2);
                dist2 = sqrt((x(i+1) - ox(i))^2 + (y(i+1) - oy(i))^2 + (z(i+1) - oz(i))^2);
                
                % 计算相邻航迹点之间的距离
                dis3 = sqrt((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2 + (z(i+1) - z(i))^2);
                
                % 惩罚项3：确保相邻两点到圆心的距离相等
                Penalty3 = Penalty1 + abs(dist1 - dist2);
                
                % 惩罚项4：确保两点之间的距离小于半径的二倍
                radius = dist1; % 半径
                if dis3 >= 2 * radius
                    Penalty4 = Penalty4 + (dis3 - 2 * radius); % 超过部分的惩罚
                end

                % 计算圆的周长
                circumference = 2 * pi * radius;
                
                % 累加周长
                totalCircumference = totalCircumference + circumference;

                radii(i) = radius;

                % 惩罚项5
                if radius < radiusMin
                    Penalty5 = Penalty5 + penaltySmall * (radiusMin - radius);
                elseif radius > radiusMax
                    Penalty5 = Penalty5 + penaltyLarge * (radius - radiusMax);
                end
            end

        PopObj(1,1)=totalCircumference; %目标值1，圆的周长之和   
        % 计算半径均值
        meanRadius = mean(radii);
        % 计算均值的倒数
            if meanRadius > 0
                reciprocalMeanRadius = 1 / meanRadius;
            else
                reciprocalMeanRadius = Inf; % 如果均值为0，倒数设为无穷大
            end
        PopObj(1,2) = reciprocalMeanRadius;%目标值1，圆的半径均值倒数 


        if PopObj(1,1)>40
            Penalty6=Penalty6+(PopObj(1,1)-40);%惩罚项6，限制总航迹长度不能大于5（500）
        end

        for i=2:nNodes                                                     % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            if delta_z/norm(a_i)>1
            Penalty7=Penalty7+delta_z/norm(a_i);%惩罚项7，限制爬升/俯冲角tan值小于1
            end
        end

        PopObj(1,1)=PopObj(1,1)+0.07*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6+Penalty7);
        PopObj(1,2)=PopObj(1,2)+0.07*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6+Penalty7);
        
        end
end
 


        %% Generate points on the Pareto front
        % function R = GetOptimum(obj,N)
        % % 
        %          R=load('SMEAGLT5.pf');
        %          % R=[20,1,0.0001];
        % % 
        % end

     
 end
end
