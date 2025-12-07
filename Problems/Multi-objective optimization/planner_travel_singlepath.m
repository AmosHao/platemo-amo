classdef planner_travel_singlepath < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            %场景设置
            dot1=[0.2,0.2,0.8];
            dot2=[2.5,6.8,0.8];
            % dot2=[3.5,3,0.8];
            D_restrain=0.4;z_min=0.8;z_max=2.4;
            
            distanse_dots=sqrt((dot1(1)-dot2(1)).^2+(dot1(2)-dot2(2)).^2+(dot1(3)-dot2(3)).^2);
            n_dots=ceil(distanse_dots/D_restrain);%向上取整
            if isempty(obj.D); obj.D =n_dots*3; end
            D=obj.D;
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
            dot1=[0.2,0.2,0.8];
            dot2=[2.5,6.8,0.8];
            % dot2=[3.5,3,0.8];
        D_min=0.4;z_min=0.8;z_max=2.4;
            %长方形
            xLimits = [2.6 4.2];
            yLimits = [0.8 2.8];
            %平行四边形
            polyX = [1.8, 1.8, 5, 5]; % x 坐标
            polyY = [3.3, 3.6, 5.3, 5]; % y 坐标
            %圆形
            circleCenter = [1.2, 3.2];
            radius = 0.6;



        distanse_dots=sqrt((dot1(1)-dot2(1)).^2+(dot1(2)-dot2(2)).^2+(dot1(3)-dot2(3)).^2);
        if size(PopDec,1)>1
        for j = 1:N 
        Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;
        xStart=dot1; xEnd=dot2;                              % 起点和终点
        z=PopDec(j,2*n+1:3*n); y=PopDec(j,n+1:2*n); x=PopDec(j,1:n);
        x=[xStart(1),x,xEnd(1)]; y=[xStart(2),y,xEnd(2)]; z=[xStart(3),z,xEnd(3)]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        nNodes=n+2;
        
        dis=zeros(1,nNodes-1);
        for i=1:nNodes-1
            if x(i+1)-x(i)<0
                Penalty1=Penalty1+(x(i)-x(i+1));%约束前进
            end
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)>D_min+0.1%比D_min高0.1 %约束最大航迹段长度
                Penalty2=Penalty2+dis(i);
            end
            if dis(i)<0.1 %约束最小航迹段长度0.1
                Penalty6=Penalty6+dis(i);
            end
        end
        PopObj(j,1)=sum(dis);   % 计算全部航迹的长度，即优化目标1
        if PopObj(j,1)>2*distanse_dots
            Penalty3=Penalty3+PopObj(j,1);
        end
        
        r0=zeros(1,nNodes-2);
        for i=1:nNodes-2
            a_i=[x(i+1)-x(i),y(i+1)-y(i)]; b_i=[x(i+2)-x(i+1),y(i+2)-y(i+1)]; r0(1,i)=dot(a_i,b_i)/(norm(a_i)*norm(b_i)+0.000001);
            if r0(1,i)<0.5
                Penalty4=(Penalty4+(0.5-r0(1,i)));%转弯角约束
            end
        end
        

        for i=2:nNodes          % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            if delta_z/norm(a_i)>1
            Penalty5=Penalty5+delta_z/norm(a_i);
            end
        end

        PopObj(j,2)=sum(z-z_min);  % 计算优化目标2  

        %计算目标三环境成本（圆心和矩心坐标设置）
        % for i=1:n
        %     dis_yuan(i)=sqrt((x(i+1)-1.2).^2+(y(i+1)-3.2).^2+(z(i+1)-0).^2);
        %     dis_juxing(i)=sqrt((x(i+1)-3.2).^2+(y(i+1)-1.8).^2+(z(i+1)-0).^2);
        % end
        % PopObj(j,3)=sum(dis_juxing)+sum(dis_yuan);
         
         dis_yuan=[];
         dis_juxing=[];
         touying=[];
         for i = 1:n
            % 当前点坐标
            currentX = x(i + 1);
            currentY = y(i + 1);
            currentZ = z(i + 1);
            % 计算点到圆的距离
            disyuan = sqrt((currentX - circleCenter(1))^2 + (currentY - circleCenter(2))^2);
            if disyuan<radius
                dis_yuan(i)=sqrt((currentX-circleCenter(1))^2+(currentY-circleCenter(2))^2+(currentZ-0)^2);
            end
            % 计算点至长方形的最近距离
            if (currentX >= xLimits(1) && currentX <= xLimits(2) && currentY >= yLimits(1) && currentY <= yLimits(2))
                distToRect = 0; % 点在长方形内部
            else
                dis_juxing(i)=sqrt((currentX-(xLimits(1)+xLimits(2))/2).^2+(currentY-(yLimits(1)+yLimits(2))/2).^2+(currentZ-0).^2);
            end
           % 计算点到平行四边形的最小距离
           % 使用 inpolygon 函数判断点是否在多边形内部
            inside = inpolygon(currentX, currentY, polyX, polyY);
            % 如果点在多边形内部，将 touying(i) 加 1
            if inside
                touying(i) = 1;
            end
            % distToParallelogram = distancePointToPolygon(currentX, currentY, polyX, polyY);    
         end

        PopObj(j,3)=1000/(sum(dis_yuan)+sum(dis_juxing)+sum(touying));

        PopObj(j,1)=PopObj(j,1)+(Penalty1+2*Penalty2+Penalty3+Penalty4+Penalty5/5+Penalty6);
        PopObj(j,2)=PopObj(j,2)+(Penalty1+2*Penalty2+Penalty3+Penalty4+Penalty5/5+Penalty6);
        PopObj(j,3)=PopObj(j,3)+(Penalty1+2*Penalty2+Penalty3+Penalty4+Penalty5/5+Penalty6);
        end
       
        else
        Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;
        xStart=dot1; xEnd=dot2;                              % 起点和终点
        z=PopDec(1,2*n+1:3*n); y=PopDec(1,n+1:2*n); x=PopDec(1,1:n);
        x=[xStart(1),x,xEnd(1)]; y=[xStart(2),y,xEnd(2)]; z=[xStart(3),z,xEnd(3)]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        nNodes=n+2;
        
        dis=zeros(1,nNodes-1);
        for i=1:nNodes-1
            if x(i+1)-x(i)<0
                Penalty1=Penalty1+(x(i)-x(i+1));%约束前进
            end
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)>D_min+0.1%比D_min高0.1 %约束最da航迹段长度
                Penalty2=Penalty2+dis(i);
            end
            if dis(i)<0.1 %约束最小航迹段长度0.1
                Penalty6=Penalty6+dis(i);
            end
        end
        PopObj(1,1)=sum(dis);   % 计算全部航迹的长度，即优化目标1
        if PopObj(1,1)>2*distanse_dots
            Penalty3=Penalty3+PopObj(1,1);
        end
        
        r0=zeros(1,nNodes-2);
        for i=1:nNodes-2
            a_i=[x(i+1)-x(i),y(i+1)-y(i)]; b_i=[x(i+2)-x(i+1),y(i+2)-y(i+1)]; r0(1,i)=dot(a_i,b_i)/(norm(a_i)*norm(b_i)+0.000001);
            if r0(1,i)<0.5
                Penalty4=Penalty4+(0.5-r0(1,i));%转弯角约束
            end
        end
        

        for i=2:nNodes          % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            if delta_z/norm(a_i)>1
            Penalty5=Penalty5+delta_z/norm(a_i);
            end
        end

        PopObj(1,2)=sum(z-z_min);  % 计算优化目标2  

        % %计算目标三环境成本（圆心和矩心坐标设置）
        % for i=1:n
        %     dis_yuan(i)=sqrt((x(i+1)-1.2).^2+(y(i+1)-3.2).^2+(z(i+1)-0).^2);
        %     dis_juxing(i)=sqrt((x(i+1)-3.2).^2+(y(i+1)-1.8).^2+(z(i+1)-0).^2);
        % end
        % PopObj(1,3)=sum(dis_juxing)+sum(dis_yuan);

         dis_yuan=[];
         dis_juxing=[];
         touying=[];
         for i = 1:n
            % 当前点坐标
            currentX = x(i + 1);
            currentY = y(i + 1);
            currentZ = z(i + 1);
            % 计算点到圆的距离
            disyuan = sqrt((currentX - circleCenter(1))^2 + (currentY - circleCenter(2))^2);
            if disyuan<radius
                dis_yuan(i)=sqrt((currentX-circleCenter(1)).^2+(currentY-circleCenter(2)).^2+(currentZ-0).^2);
            end
            % 计算点至长方形的最近距离
            if (currentX >= xLimits(1) && currentX <= xLimits(2) && currentY >= yLimits(1) && currentY <= yLimits(2))
                distToRect = 0; % 点在长方形内部
            else
                dis_juxing(i)=sqrt((currentX-(xLimits(1)+xLimits(2))/2).^2+(currentY-(yLimits(1)+yLimits(2))/2).^2+(currentZ-0).^2);
            end
           % 计算点到平行四边形的最小距离
           % 使用 inpolygon 函数判断点是否在多边形内部
            inside = inpolygon(currentX, currentY, polyX, polyY);
            % 如果点在多边形内部，将 touying(i) 加 1
            if inside
                touying(i) = 1;
            end
            % distToParallelogram = distancePointToPolygon(currentX, currentY, polyX, polyY);    
        end
        PopObj(1,3)=1000/(sum(dis_yuan)+sum(dis_juxing)+sum(touying));

        PopObj(1,1)=PopObj(1,1)+(Penalty1+2*Penalty2+Penalty3+Penalty4+Penalty5/5+Penalty6);
        PopObj(1,2)=PopObj(1,2)+(Penalty1+2*Penalty2+Penalty3+Penalty4+Penalty5/5+Penalty6);
        PopObj(1,3)=PopObj(1,3)+(Penalty1+2*Penalty2+Penalty3+Penalty4+Penalty5/5+Penalty6);

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
