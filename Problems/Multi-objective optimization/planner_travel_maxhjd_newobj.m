classdef planner_travel_maxhjd_newobj < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            % if isempty(obj.D); obj.D = 600; end
            %场景设置
            if isempty(obj.dots); obj.dots = [7000,3000,800,9500,6000,800]; end
            % if isempty(obj.dots); obj.dots = [7000,3000,800,9500,6000,800]; end
            dots=obj.dots;
            dot1=[dots(1,1),dots(1,2),dots(1,3)];
            dot2=[dots(1,4),dots(1,5),dots(1,6)];
            % dot1=[500,500,800];
            % % dot2=[1800,8950,800];
            % dot2=[7000,3000,800];
            D_restrain=100;z_min=800;z_max=2400;
            % n_dots=obj.D/3;
            disx=abs(dot2(1,1)-dot1(1,1));
            disy=abs(dot2(1,2)-dot1(1,2));
            dis_min=min(disx,disy);
            distanse_dots=sqrt((dot1(1)-dot2(1)).^2+(dot1(2)-dot2(2)).^2+(dot1(3)-dot2(3)).^2);
            % n_dots=ceil(1.1*dis_min/D_restrain);%向上取整
            n_dots=ceil(1.5*distanse_dots/D_restrain);%向上取整
            if isempty(obj.D); obj.D =n_dots*3; end
            % D=obj.D;
            obj.lower(1:n_dots) = min(dot1(1),dot2(1))*ones(1,n_dots);
            obj.upper(1:n_dots) = max(dot1(1),dot2(1))*ones(1,n_dots);
            obj.lower(n_dots+1:2*n_dots) = min(dot1(2),dot2(2))*ones(1,n_dots);
            obj.upper(n_dots+1:2*n_dots) = max(dot1(2),dot2(2))*ones(1,n_dots);
            obj.lower(2*n_dots+1:3*n_dots) = z_min*ones(1,n_dots);
            obj.upper(2*n_dots+1:3*n_dots) = z_max*ones(1,n_dots);

            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        M=obj.M;
        n=D/3;
        PopObj =zeros(N,M);
        %场景设置
            % dot1=[500,500,800];
            % dot2=[7000,3000,8000];
                        dots=obj.dots;
            dot1=[dots(1,1),dots(1,2),dots(1,3)];
            dot2=[dots(1,4),dots(1,5),dots(1,6)];
            D_restrain=100;z_min=800;z_max=2400;
            %长方形
            xLimits = [6500 ,8000;9100 ,10000;2000, 3100;6800, 8100;1900,3300];
            yLimits = [7000 ,8700;6300 ,8600;1200, 2200;3300, 4400;6300,8300];
            %平行四边形
            polyX = [0, 0, 2000, 2000]; % x 坐标
            polyY = [7300, 7450, 6050,5900]; % y 坐标
            %圆形
            circleCenter = [4600, 4700;4500,2200;4800,8700];
            radius = [900;600;600];
        distanse_dots=sqrt((dot1(1)-dot2(1)).^2+(dot1(2)-dot2(2)).^2+(dot1(3)-dot2(3)).^2);
        if size(PopDec,1)>1
        for j = 1:N 
        Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;
        % xStart=dot1; xEnd=dot2;                              % 起点和终点
        z=PopDec(j,2*n+1:3*n); y=PopDec(j,n+1:2*n); x=PopDec(j,1:n);
        % x=[min(dot1(1),dot2(1)),x,max(dot1(1),dot2(1))]; y=[min(dot1(2),dot2(2)),y,max(dot1(2),dot2(2))]; z=[min(dot1(3),dot2(3)),z,max(dot1(3),dot2(3))]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        x=[dot1(1),x,dot2(1)]; y=[dot1(2),y,dot2(2)]; z=[dot1(3),z,dot2(3)]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        nNodes=n+2;
        
        dis=zeros(1,nNodes-1);
        for i=1:nNodes-1
            if dot1(1,1)<dot2(1,1)
                if x(i+1)-x(i)<0
                    Penalty1=Penalty1+(x(i)-x(i+1));%约束前进
                end
            else
               if x(i+1)-x(i)>0
                Penalty1=Penalty1+(x(i+1)-x(i));%约束前进
               end 
            end
            if dot1(1,2)<dot2(1,2)
            if y(i+1)-y(i)<0
                Penalty1=Penalty1+(y(i)-y(i+1));%约束前进
            end
            else
            if y(i+1)-y(i)>0
                Penalty1=Penalty1+(y(i+1)-y(i));%约束前进
            end   
            end
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)>D_restrain%最大航迹段长度
                Penalty2=Penalty2+dis(i);
            end
            if dis(i)<5 %约束最小航迹段长度5
                Penalty6=Penalty6+dis(i);
            end
        end
        PopObj(j,1)=sum(dis);   % 计算全部航迹的长度，即优化目标1
        if PopObj(j,1)>1.5*distanse_dots
            Penalty3=Penalty3+PopObj(j,1);
        end
        
        r0=zeros(1,nNodes-2);
        for i=1:nNodes-2
            a_i=[x(i+1)-x(i),y(i+1)-y(i)]; b_i=[x(i+2)-x(i+1),y(i+2)-y(i+1)]; 
             r0(1,i)=dot(a_i,b_i)/(norm(a_i)*norm(b_i)+0.000001);
            % r0(1, i) = abs(dot(a_i, b_i)) / (norm(a_i) * norm(b_i) + 0.000001);
            if r0(1,i)<0.5
                Penalty4=(Penalty4+(0.5-r0(1,i)));%转弯角约束
            end
        end
        
        f0=zeros(1,nNodes-1);
        for i=2:nNodes          % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            f0(1,i)=delta_z/(norm(a_i)+0.000001);
            if delta_z/norm(a_i)>1
            Penalty5=Penalty5+delta_z/(norm(a_i)+0.000001);
            end
        end

        % PopObj(j,2)=sum(z-z_min);  % 计算优化目标2  
        r1=ones(1,nNodes-2);
        r2=r1-r0;
        f1=zeros(1,nNodes-2);
        for i=1:nNodes-2
            f1(1,i)=abs(f0(1,i+1)-f0(1,i));
        end
        % a=sum(r2);
        PopObj(j,2)=sum(r2)+sum(f1);

        %计算目标三环境成本（圆心和矩心坐标设置）
      r_ju=[];
      for i=1:size(xLimits,1)
          r_ju(i,:)=sqrt((xLimits(i,1) - xLimits(i,2))^2 + (yLimits(i,1) - yLimits(i,2))^2)/2;
      end
         
         dis_yuan=[];
         dis_juxing=[];
         obj_yuan=[];
         obj_juxing=[];
         for i=1:nNodes-2
             for k=1:size(xLimits,1)
             dis_juxing(k,i)=sqrt(((xLimits(k,1)-xLimits(k,2))/2-x(i+1))^2 +((yLimits(k,1)-yLimits(k,2))/2-y(i+1))^2);  
             obj_juxing(k,i)=exp(-dis_juxing(k,i)^2/r_ju(k,1)^2);
             end
         end
         for i=1:nNodes-2
             for k=1:size(circleCenter,1)
             dis_yuan(k,i)=sqrt((circleCenter(k,1)-x(i+1))^2 +(circleCenter(k,2)-y(i+1))^2);  
             obj_yuan(k,i)=exp(-dis_yuan(k,i)^2/radius(k,1)^2);
             end
         end
        PopObj(j,3)=sum(sum(obj_yuan)+sum(obj_juxing));
Penalty4=5000*Penalty4;
Penalty5=5000*Penalty5;
PopObj(j,1)=PopObj(j,1)/1000;
        PopObj(j,1)=PopObj(j,1)+(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(j,2)=PopObj(j,2)+(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(j,3)=PopObj(j,3)+(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        end
       
        else
            j=1;
        Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;
        % xStart=dot1; xEnd=dot2;                              % 起点和终点
        z=PopDec(j,2*n+1:3*n); y=PopDec(j,n+1:2*n); x=PopDec(j,1:n);
        % x=[min(dot1(1),dot2(1)),x,max(dot1(1),dot2(1))]; y=[min(dot1(2),dot2(2)),y,max(dot1(2),dot2(2))]; z=[min(dot1(3),dot2(3)),z,max(dot1(3),dot2(3))]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        x=[dot1(1),x,dot2(1)]; y=[dot1(2),y,dot2(2)]; z=[dot1(3),z,dot2(3)]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        nNodes=n+2;
        
        dis=zeros(1,nNodes-1);
        for i=1:nNodes-1
            if dot1(1,1)<dot2(1,1)
                if x(i+1)-x(i)<0
                    Penalty1=Penalty1+(x(i)-x(i+1));%约束前进
                end
            else
               if x(i+1)-x(i)>0
                Penalty1=Penalty1+(x(i+1)-x(i));%约束前进
               end 
            end
            if dot1(1,2)<dot2(1,2)
            if y(i+1)-y(i)<0
                Penalty1=Penalty1+(y(i)-y(i+1));%约束前进
            end
            else
            if y(i+1)-y(i)>0
                Penalty1=Penalty1+(y(i+1)-y(i));%约束前进
            end   
            end
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)>D_restrain%最大航迹段长度
                Penalty2=Penalty2+dis(i);
            end
            if dis(i)<5 %约束最小航迹段长度5
                Penalty6=Penalty6+dis(i);
            end
        end
        PopObj(j,1)=sum(dis);   % 计算全部航迹的长度，即优化目标1
        if PopObj(j,1)>1.5*distanse_dots
            Penalty3=Penalty3+PopObj(j,1);
        end
        
        r0=zeros(1,nNodes-2);
        for i=1:nNodes-2
            a_i=[x(i+1)-x(i),y(i+1)-y(i)]; b_i=[x(i+2)-x(i+1),y(i+2)-y(i+1)]; 
             r0(1,i)=dot(a_i,b_i)/(norm(a_i)*norm(b_i)+0.000001);
            % r0(1, i) = abs(dot(a_i, b_i)) / (norm(a_i) * norm(b_i) + 0.000001);
            if r0(1,i)<0.5
                Penalty4=(Penalty4+(0.5-r0(1,i)));%转弯角约束
            end
        end
        
        f0=zeros(1,nNodes-1);
        for i=2:nNodes          % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            f0(1,i)=delta_z/(norm(a_i)+0.000001);
            if delta_z/norm(a_i)>1
            Penalty5=Penalty5+delta_z/(norm(a_i)+0.000001);
            end
        end

        % PopObj(j,2)=sum(z-z_min);  % 计算优化目标2  
        r1=ones(1,nNodes-2);
        r2=r1-r0;
        f1=zeros(1,nNodes-2);
        for i=1:nNodes-2
            f1(1,i)=abs(f0(1,i+1)-f0(1,i));
        end
        % a=sum(r2);
        PopObj(j,2)=sum(r2)+sum(f1);

        %计算目标三环境成本（圆心和矩心坐标设置）
      r_ju=[];
      for i=1:size(xLimits,1)
          r_ju(i,:)=sqrt((xLimits(i,1) - xLimits(i,2))^2 + (yLimits(i,1) - yLimits(i,2))^2)/2;
      end
         
         dis_yuan=[];
         dis_juxing=[];
         obj_yuan=[];
         obj_juxing=[];
         for i=1:nNodes-2
             for k=1:size(xLimits,1)
             dis_juxing(k,i)=sqrt(((xLimits(k,1)-xLimits(k,2))/2-x(i+1))^2 +((yLimits(k,1)-yLimits(k,2))/2-y(i+1))^2);  
             obj_juxing(k,i)=exp(-dis_juxing(k,i)^2/r_ju(k,1)^2);
             end
         end
         for i=1:nNodes-2
             for k=1:size(circleCenter,1)
             dis_yuan(k,i)=sqrt((circleCenter(k,1)-x(i+1))^2 +(circleCenter(k,2)-y(i+1))^2);  
             obj_yuan(k,i)=exp(-dis_yuan(k,i)^2/radius(k,1)^2);
             end
         end
        PopObj(j,3)=sum(sum(obj_yuan)+sum(obj_juxing));
Penalty4=5000*Penalty4;
Penalty5=5000*Penalty5;
PopObj(j,1)=PopObj(j,1)/1000;
        PopObj(j,1)=PopObj(j,1)+(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(j,2)=PopObj(j,2)+(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(j,3)=PopObj(j,3)+(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);

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
 function dist = distancePointToPolygon(px, py, polyX, polyY)
    % polyX, polyY 是平行四边形的顶点坐标
    % 使用 MATLAB 的内置方法来进行多边形测试
    % 计算到每个边的距离
    distances = [];
    for j = 1:length(polyX)
        % 得到当前边的两个顶点
        pt1 = [polyX(j), polyY(j)];
        pt2 = [polyX(mod(j, length(polyX))+1), polyY(mod(j, length(polyY))+1)];
        
        % 计算边的分量
        edge = pt2 - pt1;
        lengthEdge = norm(edge);
        
        if lengthEdge == 0
            % 如果边的长度为 0，使用点到顶点的距离
            distances(end + 1) = norm([px, py] - pt1);
        else
            % 计算当前点到边的投影
            projection = (dot([px, py] - pt1, edge) / lengthEdge^2) * edge + pt1;       
            % 计算距离
            distances(end + 1) = norm([px, py] - projection);
        end
    end
    
    % 返回最小距离
    dist = min(distances);
end
