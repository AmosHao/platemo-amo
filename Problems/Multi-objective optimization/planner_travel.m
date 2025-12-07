classdef planner_travel < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D =59; end
            D=obj.D;
            % n=D/3; 
            obj.lower(1:5)=1*ones(1,5);
            obj.upper(1:5)=5*ones(1,5);
            obj.encoding(1:5) = 5*ones(1,5);

            obj.lower(6:23) = zeros(1,18);
            obj.upper(6:23) = 5*ones(1,18);
            obj.lower(24:41) = zeros(1,18);
            obj.upper(24:41) = 7*ones(1,18);
            obj.lower(42:59) = 0.8*ones(1,18);
            obj.upper(42:59) = 2.4*ones(1,18);
            obj.encoding(6:59) = ones(1,54);
         
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        M=obj.M;
        PopObj =zeros(N,M);
        % 坐标映射表 (坐标点与编号的关系)
coordinates = [2.5, 0.5, 0.8, 0.2; 
               3.5, 3, 0.8, 0.5; 
               2.5, 5, 0.8, 0.8; 
               2.5, 6.8, 0.8, 0.3; 
               0.5, 5, 0.8, 1.0];  % 5个客户点的坐标和货物重量
if size(PopDec,1)>1
     for j = 1:N 
        Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0; Penalty6=0;
        xStart=[0.2,0.2,0.8];   % 起点
        indices = PopDec(j, 1:5); % 前五列的编号
        % 提取每个编号对应的坐标
        coords = zeros(5, 3); % 5个点，每个点有3个坐标
        for k = 1:5
            coords(k, :) = coordinates(indices(k), 1:3);
        end
        %将客户点坐标插入到航迹点坐标中
        z=PopDec(j,42:59); y=PopDec(j,24:41); x=PopDec(j,6:23);
        x=[xStart(1),x(1:3),coords(1,1),x(4:6),coords(2,1),x(7:9),coords(3,1),x(10:12),coords(4,1),x(13:15),coords(5,1),x(16:18),xStart(1)];
        y=[xStart(2),y(1:3),coords(1,2),y(4:6),coords(2,2),y(7:9),coords(3,2),y(10:12),coords(4,2),y(13:15),coords(5,2),y(16:18),xStart(2)];
        z=[xStart(3),z(1:3),coords(1,3),z(4:6),coords(2,3),z(7:9),coords(3,3),z(10:12),coords(4,3),z(13:15),coords(5,3),z(16:18),xStart(3)];
        %放缩坐标轴
        x=x/2;
        y=y/2;
        z=z/0.15;
        %惩罚项6，使每段航迹点靠近客户
        d1=0;
        for i=2:4
        d1=d1+sqrt((x(5)-x(i)).^2+(y(5)-y(i)).^2+(z(5)-z(i)).^2);
        end
        d2=0;
        for i=6:8
        d2=d2+sqrt((x(9)-x(i)).^2+(y(9)-y(i)).^2+(z(9)-z(i)).^2);
        end
        d3=0;
        for i=10:12
        d3=d3+sqrt((x(13)-x(i)).^2+(y(13)-y(i)).^2+(z(13)-z(i)).^2);
        end
        d4=0;
        for i=14:16
        d4=d4+sqrt((x(17)-x(i)).^2+(y(17)-y(i)).^2+(z(17)-z(i)).^2);
        end        
        d5=0;
        for i=18:20
        d5=d5+sqrt((x(21)-x(i)).^2+(y(21)-y(i)).^2+(z(21)-z(i)).^2);
        end
        d6=0;
        for i=22:24
        d6=d6+sqrt((x(25)-x(i)).^2+(y(25)-y(i)).^2+(z(25)-z(i)).^2);
        end
        dd=d1+d2+d3+d4+d5+d6;
        Penalty6=dd;

        nleg=25;%共25个航段
        dis=zeros(1,nleg);
        for i=1:nleg-1
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)<0.1
                Penalty1=Penalty1+dis(i);%惩罚项1，限制两点之间距离大于0.1
            end
        end
        PopObj(j,1)=sum(dis); %目标值1，总航迹长度                                              % 计算全部航迹的长度，即优化目标1
        if PopObj(j,1)>10
            Penalty2=Penalty2+PopObj(j,1);%惩罚项2，限制总航迹长度不能大于10
        end
        r0=zeros(1,nleg-2);
        for i=1:nleg-2
            a_i=[x(i+1)-x(i),y(i+1)-y(i)]; b_i=[x(i+2)-x(i+1),y(i+2)-y(i+1)]; r0(1,i)=dot(a_i,b_i)/(norm(a_i)*norm(b_i)+0.000001);
            if r0(1,i)<0.5
                Penalty3=Penalty3+(0.5-r0(1,i));%惩罚项3，限制转弯角余弦值大于0.5
            end
        end
        
        for i=2:nleg  % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            if delta_z/norm(a_i)>1
            Penalty4=Penalty4+delta_z/norm(a_i);%惩罚项4，限制爬升/俯冲角tan值小于1
            end
        end
        %速度约束km/h
        v_m=83; 
        v_up_m=21.6;
        v_down_m=14.4;
        v_metrics=zeros(2,6);%存储6个大航段的水平和上升速度
        % 提取每个编号对应的货物重量
        coords_weight = zeros(5, 1); 
        for k =1:5
          coords_weight(k, :) = coordinates(indices(k), 4);
        end
        for i=1:5
          v_metrics(1,i)=v_m-3*(sum(coords_weight(i:5, :)));
          v_metrics(2,i)=v_up_m-2*(sum(coords_weight(i:5, :)));
        end
        v_metrics(1,6)=v_m;
        v_metrics(2,6)=v_up_m;

       dis1=sum(dis(1:5));
       dis2=sum(dis(5:9));
       dis3=sum(dis(9:13));
       dis4=sum(dis(13:17));
       dis5=sum(dis(17:21));
       dis6=sum(dis(21:25));
        t1 = dis1/v_metrics(1,1)+dis2/v_metrics(1,2)+dis3/v_metrics(1,3)+dis4/v_metrics(1,4)+dis5/v_metrics(1,5)+dis6/v_metrics(1,6);
        t2 = sum(0.04 ./ v_metrics(2, :));
        t3=6*0.04/v_down_m;
        PopObj(j,2)=t1+t2+t3;%目标函数2，时间
        %目标函数3
        for i=1:nleg-1
            dis_yuan(i)=sqrt((x(i)-1.2).^2+(y(i)-3.2).^2+(z(i)-0).^2);
            dis_juxing(i)=sqrt((x(i)-3.2).^2+(y(i)-1.8).^2+(z(i)-0).^2);
        end
        PopObj(j,3)=sum(dis_juxing)+sum(dis_yuan);

        % 期望的数字集合
        validNumbers = [1, 2, 3, 4, 5];
        % 提取前五列
        currentValues = PopDec(j, 1:5);
        %惩罚项5，必须遍历每个客户点
        if ~isequal(sort(currentValues), validNumbers)
           Penalty5=sum(sort(currentValues)-validNumbers).^2;
        end
        PopObj(j,1)=PopObj(j,1)+8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(j,2)=PopObj(j,2)+8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(j,3)=PopObj(j,3)+8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
    end
% PopObj(i,1)=sum(sqrt((PopDec(i,2:n)-PopDec(i,1:n-1)).^2+(PopDec(i,n+2:2*n)-PopDec(i,n+1:2*n-1)).^2+(PopDec(i,2*n+2:3*n)-PopDec(i,2*n+1:3*n-1)).^2+(PopDec(i,1)-p1).^2+(PopDec(i,n+1)-p2).^2+(PopDec(i,2*n+1)-p3).^2+(PopDec(i,n)-q1).^2+(PopDec(i,2*n)-q2).^2+(PopDec(i,3*n)-q3).^2));
% g1=1-(((PopDec(i,1)-p1)*(PopDec(i,2)-PopDec(i,1)))+((PopDec(i,n+1)-p2)*(PopDec(i,n+2)-PopDec(i,n+1))))/sqrt(((PopDec(i,1)-p1).^2+(PopDec(i,n+1)-p2).^2)*((PopDec(i,2)-PopDec(i,1)).^2+(PopDec(i,n+2)-PopDec(i,n+1)).^2));
% g2=1-(((PopDec(i,n)-PopDec(i,n-1))*(q1-PopDec(i,n)))+((PopDec(i,2*n)-PopDec(i,2*n-1))*(q2-PopDec(i,2*n))))/sqrt(((PopDec(i,n)-PopDec(i,n-1)).^2+(PopDec(i,2*n)-PopDec(i,2*n-1)).^2)*((q1-PopDec(i,n)).^2+(q2-PopDec(i,2*n)).^2));
% PopObj(i,2)=sum(PopDec(i,2*n+1:3*n));
% 
% PopObj(i,3)=g1+g2+sum(1-(((PopDec(i,2:n-1)-PopDec(i,1:n-2)).*(PopDec(i,3:n)-PopDec(i,2:n-1)))+((PopDec(i,n+2:2*n-1)-PopDec(i,n+1:2*n-2)).*(PopDec(i,n+3:2*n)-PopDec(i,n+2:2*n-1))))/sqrt(((PopDec(i,2:n-1)-PopDec(i,1:n-2)).^2+(PopDec(i,n+2:2*n-1)-PopDec(i,n+1:2*n-2)).^2).*((PopDec(i,3:n)-PopDec(i,2:n-1)).^2+(PopDec(i,n+3:2*n)-PopDec(i,n+2:2*n-1)).^2)));
   
            
else
        Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;
        xStart=[0.2,0.2,0.8];   % 起点
        indices = PopDec(1, 1:5); % 前五列的编号
        % 提取每个编号对应的坐标
        coords = zeros(5, 3); % 5个点，每个点有3个坐标
        for k = 1:5
            coords(k, :) = coordinates(indices(k), 3);
        end
        %将客户点坐标插入到航迹点坐标中
        z=PopDec(1,42:59); y=PopDec(1,24:41); x=PopDec(1,6:23);
        x=[xStart(1),x(1:3),coords(1,1),x(4:6),coords(2,1),x(7:9),coords(3,1),x(10:12),coords(4,1),x(13:15),coords(5,1),x(16:18),xStart(1)];
        y=[xStart(2),y(1:3),coords(1,2),y(4:6),coords(2,2),y(7:9),coords(3,2),y(10:12),coords(4,2),y(13:15),coords(5,2),y(16:18),xStart(2)];
        z=[xStart(3),z(1:3),coords(1,3),z(4:6),coords(2,3),z(7:9),coords(3,3),z(10:12),coords(4,3),z(13:15),coords(5,3),z(16:18),xStart(3)];
        %放缩坐标轴
        x=x/2;
        y=y/2;
        z=z/0.15;

        %惩罚项6，使每段航迹点靠近客户
        d1=0;
        for i=2:4
        d1=d1+sqrt((x(5)-x(i)).^2+(y(5)-y(i)).^2+(z(5)-z(i)).^2);
        end
        d2=0;
        for i=6:8
        d2=d2+sqrt((x(9)-x(i)).^2+(y(9)-y(i)).^2+(z(9)-z(i)).^2);
        end
        d3=0;
        for i=10:12
        d3=d3+sqrt((x(13)-x(i)).^2+(y(13)-y(i)).^2+(z(13)-z(i)).^2);
        end
        d4=0;
        for i=14:16
        d4=d4+sqrt((x(17)-x(i)).^2+(y(17)-y(i)).^2+(z(17)-z(i)).^2);
        end        
        d5=0;
        for i=18:20
        d5=d5+sqrt((x(21)-x(i)).^2+(y(21)-y(i)).^2+(z(21)-z(i)).^2);
        end
        d6=0;
        for i=22:24
        d6=d6+sqrt((x(25)-x(i)).^2+(y(25)-y(i)).^2+(z(25)-z(i)).^2);
        end
        dd=d1+d2+d3+d4+d5+d6;
        Penalty6=dd*10;

        nleg=25;%共25个航段
        dis=zeros(1,nleg);
        for i=1:nleg-1
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)<0.1
                Penalty1=Penalty1+dis(i);%惩罚项1，限制两点之间距离大于0.1
            end
        end
        PopObj(1,1)=sum(dis); %目标值1，总航迹长度                                              % 计算全部航迹的长度，即优化目标1
        if PopObj(1,1)>10
            Penalty2=Penalty2+PopObj(1,1);%惩罚项2，限制总航迹长度不能大于10
        end
        r0=zeros(1,nleg-2);
        for i=1:nleg-2
            a_i=[x(i+1)-x(i),y(i+1)-y(i)]; b_i=[x(i+2)-x(i+1),y(i+2)-y(i+1)]; r0(1,i)=dot(a_i,b_i)/(norm(a_i)*norm(b_i)+0.000001);
            if r0(1,i)<0.5
                Penalty3=Penalty3+(0.5-r0(1,i));%惩罚项3，限制转弯角余弦值大于0.5
            end
        end
        
        for i=2:nleg  % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            if delta_z/norm(a_i)>1
            Penalty4=Penalty4+delta_z/norm(a_i);%惩罚项4，限制爬升/俯冲角tan值小于1
            end
        end
        %速度约束km/h
        v_m=83; 
        v_up_m=21.6;
        v_down_m=14.4;
        v_metrics=zeros(2,6);%存储6个大航段的水平和上升速度
        % 提取每个编号对应的货物重量
        coords_weight = zeros(5, 1); 
        for k =1:5
          coords_weight(k, :) = coordinates(indices(k), 4);
        end
        for i=1:5
          v_metrics(1,i)=v_m-3*(sum(coords_weight(i:5, :)));
          v_metrics(2,i)=v_up_m-2*(sum(coords_weight(i:5, :)));
        end
        v_metrics(1,6)=v_m;
        v_metrics(2,6)=v_up_m;

       dis1=sum(dis(1:5));
       dis2=sum(dis(5:9));
       dis3=sum(dis(9:13));
       dis4=sum(dis(13:17));
       dis5=sum(dis(17:21));
       dis6=sum(dis(21:25));
        t1 = dis1/v_metrics(1,1)+dis2/v_metrics(1,2)+dis3/v_metrics(1,3)+dis4/v_metrics(1,4)+dis5/v_metrics(1,5)+dis6/v_metrics(1,6);
        t2 = sum(0.04 ./ v_metrics(2, :));
        t3=6*0.04/v_down_m;
        PopObj(1,2)=t1+t2+t3;%目标函数2，时间
        %目标函数3
        for i=1:nleg
            dis_yuan(i)=sqrt((x(i)-1.2).^2+(y(i)-3.2).^2+(z(i)-0).^2);
            dis_juxing(i)=sqrt((x(i)-3.2).^2+(y(i)-1.8).^2+(z(i)-0).^2);
        end
        PopObj(1,3)=sum(dis_juxing)+sum(dis_yuan);

        % 期望的数字集合
        validNumbers = [1, 2, 3, 4, 5];
        % 提取前五列
        currentValues = PopDec(1, 1:5);
        %惩罚项5，必须遍历每个客户点
        if ~isequal(sort(currentValues), validNumbers)
           Penalty5=sum(sort(currentValues)-validNumbers).^2;
        end

        PopObj(1,1)=PopObj(1,1)+8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(1,2)=PopObj(1,2)+8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(1,3)=PopObj(1,3)+8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
       
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
function z=zCalculation(x,y)                                               % 利用x和y坐标计算z坐标
a=3/2*pi; b=0.1; c=0.9; d=0.5; e=0.5; f=0.3;
xb=x/15;yb=y/18;
z1=sin(0.2*yb+a)+b*sin(xb)+c*cos(d*sqrt((xb.^2+yb.^2)/5))+e*sin(e*sqrt((xb.^2+yb.^2)/5))+f*cos(yb);  % 丘陵地形的计算
x0=[50.0,100.0,100.0,130.0,160.0,70]; y0 =[60.0,160.0,100.0,20.0,100.0,30]; % 凸起山峰的中心x,y坐标，定5是为了有5个山峰
xt=[140.0,280.0,150.0,160.0,170.0,170]; yt =[20.0,220.0,280.0,190.0,230.0,150]; % 地形轮廓参数
h=[0.7,2.50,3.2,2.34,1.77,1.8];
z2=0;
for i=1:6
    z2=z2+h(i)*exp(-(x-x0(i))^2/xt(i)-(y-y0(i))^2/yt(i));                  % 山峰对应点的地形高度
end
z=max(z1,z2);                                                              % 取对应点的实际高度值
end


