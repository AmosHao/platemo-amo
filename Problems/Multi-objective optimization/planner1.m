classdef planner1 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
 methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D =30; end
            D=obj.D;
            n=D/3; 
            obj.lower(1:3*n) = zeros(1,3*n);
            obj.upper(1:2*n) = 200*ones(1,2*n);
            obj.upper(2*n+1:3*n) = 3*ones(1,n);
            obj.encoding = ones(1,obj.D);
         
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
        N=size(PopDec,1);
        D=obj.D;
        M=obj.M;
        n=D/3;
        PopObj =zeros(N,M);
        if size(PopDec,1)>1
        for j = 1:N 
                 Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;
                  xStart=[0,0,0.25]; xEnd=[200,200,3];                              % 起点和终点
        z=PopDec(j,2*n+1:3*n); y=PopDec(j,n+1:2*n); x=PopDec(j,1:n);
        x=[xStart(1),x,xEnd(1)]; y=[xStart(2),y,xEnd(2)]; z=[xStart(3),z,xEnd(3)]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        nNodes=n+2;
        
        dis=zeros(1,nNodes-1);
        for i=1:nNodes-1
            if x(i+1)-x(i)<0
                Penalty1=Penalty1+(x(i)-x(i+1));%惩罚项1，保证前进
            end
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)<2
                Penalty2=Penalty2+dis(i);%惩罚项2，限制两点之间距离大于2
            end
        end
        PopObj(j,1)=sum(dis); %目标值1，总航迹长度                                              % 计算全部航迹的长度，即优化目标1
        if PopObj(j,1)>500
            Penalty6=Penalty6+PopObj(j,1);%惩罚项6，限制总航迹长度不能大于500
        end
        r0=zeros(1,nNodes-2);
        for i=1:nNodes-2
            a_i=[x(i+1)-x(i),y(i+1)-y(i)]; b_i=[x(i+2)-x(i+1),y(i+2)-y(i+1)]; r0(1,i)=dot(a_i,b_i)/(norm(a_i)*norm(b_i)+0.000001);
            if r0(1,i)<0.5
                Penalty3=Penalty3+(0.5-r0(1,i));%惩罚项3，限制转弯角余弦值大于0.5
            end
        end
        PopObj(j,2)=sum(1-r0);                                              % 计算优化目标2  
        PopObj(j,3)=sum(z);                                                 % 优化目标3，让飞行高度尽量低

                %高度惩罚1
        HConstraint=zeros(1,10*(nNodes-1));
        height=zeros(1,10*(nNodes-1));                                     % 存储每一个航迹点（包括两个航迹点中间的采样点）的地形高度,没有管终点的情况
        for k=1:nNodes-1
            for i=0:9
                aux_x=x(k)+0.1*i*(x(k+1)-x(k)); aux_y=y(k)+0.1*i*(y(k+1)-y(k)); % 航迹点及其临时采样点的x,y坐标
               if k==1&&i==0
                   height(1,1+i+10*(k-1))=0.2;%最低飞行高度
               else
                   height(1,1+i+10*(k-1))=zCalculation(aux_x,aux_y)+0.05;       % 计算中间划分点对应的地形高度
               end
                HConstraint(1,1+i+10*(k-1))=height(1,1+i+10*(k-1))-(z(k)+0.1*i*(z(k+1)-z(k)));
                if HConstraint(1,1+i+10*(k-1))>0
                    Penalty4=Penalty4+HConstraint(1,1+i+10*(k-1));
                end
            end
        end
        
        for i=2:nNodes                                                     % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            if delta_z/norm(a_i)>1
            Penalty5=Penalty5+delta_z/norm(a_i);%惩罚项5，限制爬升/俯冲角tan值小于1
            end
        end

        PopObj(j,1)=PopObj(j,1)+10^8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(j,2)=PopObj(j,2)+10^8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(j,3)=PopObj(j,3)+10^8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
    end
% PopObj(i,1)=sum(sqrt((PopDec(i,2:n)-PopDec(i,1:n-1)).^2+(PopDec(i,n+2:2*n)-PopDec(i,n+1:2*n-1)).^2+(PopDec(i,2*n+2:3*n)-PopDec(i,2*n+1:3*n-1)).^2+(PopDec(i,1)-p1).^2+(PopDec(i,n+1)-p2).^2+(PopDec(i,2*n+1)-p3).^2+(PopDec(i,n)-q1).^2+(PopDec(i,2*n)-q2).^2+(PopDec(i,3*n)-q3).^2));
% g1=1-(((PopDec(i,1)-p1)*(PopDec(i,2)-PopDec(i,1)))+((PopDec(i,n+1)-p2)*(PopDec(i,n+2)-PopDec(i,n+1))))/sqrt(((PopDec(i,1)-p1).^2+(PopDec(i,n+1)-p2).^2)*((PopDec(i,2)-PopDec(i,1)).^2+(PopDec(i,n+2)-PopDec(i,n+1)).^2));
% g2=1-(((PopDec(i,n)-PopDec(i,n-1))*(q1-PopDec(i,n)))+((PopDec(i,2*n)-PopDec(i,2*n-1))*(q2-PopDec(i,2*n))))/sqrt(((PopDec(i,n)-PopDec(i,n-1)).^2+(PopDec(i,2*n)-PopDec(i,2*n-1)).^2)*((q1-PopDec(i,n)).^2+(q2-PopDec(i,2*n)).^2));
% PopObj(i,2)=sum(PopDec(i,2*n+1:3*n));
% 
% PopObj(i,3)=g1+g2+sum(1-(((PopDec(i,2:n-1)-PopDec(i,1:n-2)).*(PopDec(i,3:n)-PopDec(i,2:n-1)))+((PopDec(i,n+2:2*n-1)-PopDec(i,n+1:2*n-2)).*(PopDec(i,n+3:2*n)-PopDec(i,n+2:2*n-1))))/sqrt(((PopDec(i,2:n-1)-PopDec(i,1:n-2)).^2+(PopDec(i,n+2:2*n-1)-PopDec(i,n+1:2*n-2)).^2).*((PopDec(i,3:n)-PopDec(i,2:n-1)).^2+(PopDec(i,n+3:2*n)-PopDec(i,n+2:2*n-1)).^2)));
   
            
         else
                Penalty1=0; Penalty2=0; Penalty3=0; Penalty4=0; Penalty5=0;Penalty6=0;
                 xStart=[0,0,0.25]; xEnd=[200,200,3];                              % 起点和终点
        z=PopDec(1,2*n+1:3*n); y=PopDec(1,n+1:2*n); x=PopDec(1,1:n);
        x=[xStart(1),x,xEnd(1)]; y=[xStart(2),y,xEnd(2)]; z=[xStart(3),z,xEnd(3)]; % 经过这样的处理[x(1),y(1),z(1)]就是航迹的起始点
        nNodes=n+2;
        
        dis=zeros(1,nNodes-1);
        for i=1:nNodes-1
            if x(i+1)-x(i)<0
                Penalty1=Penalty1+(x(i)-x(i+1));
            end
            dis(i)=sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2+(z(i+1)-z(i)).^2); % 中间每相邻两个航迹点之间的距离
            if dis(i)<2
                Penalty2=Penalty2+dis(i);
            end
        end
        PopObj(1,1)=sum(dis);                                               % 计算全部航迹的长度，即优化目标1
        if PopObj(1,1)>500
            Penalty6=Penalty6+PopObj(1,1);
        end
        r0=zeros(1,nNodes-2);
        for i=1:nNodes-2
            a_i=[x(i+1)-x(i),y(i+1)-y(i)]; b_i=[x(i+2)-x(i+1),y(i+2)-y(i+1)]; r0(1,i)=dot(a_i,b_i)/(norm(a_i)*norm(b_i)+0.000001);
            if r0(1,i)<0.5
                Penalty3=Penalty3+(0.5-r0(1,i));
            end
        end
        PopObj(1,2)=sum(1-r0);                                              % 计算优化目标2  
        PopObj(1,3)=sum(z);                                                 % 让飞行高度尽量低
                %高度惩罚1
        HConstraint=zeros(1,10*(nNodes-1));
        height=zeros(1,10*(nNodes-1));                                     % 存储每一个航迹点（包括两个航迹点中间的采样点）的地形高度,没有管终点的情况
        for k=1:nNodes-1
            for i=0:9
                aux_x=x(k)+0.1*i*(x(k+1)-x(k)); aux_y=y(k)+0.1*i*(y(k+1)-y(k)); % 航迹点及其临时采样点的x,y坐标
               if k==1&&i==0
                   height(1,1+i+10*(k-1))=0.2;
               else
                   height(1,1+i+10*(k-1))=zCalculation(aux_x,aux_y)+0.05;       % 计算中间划分点对应的地形高度
               end
                HConstraint(1,1+i+10*(k-1))=height(1,1+i+10*(k-1))-(z(k)+0.1*i*(z(k+1)-z(k)));
                if HConstraint(1,1+i+10*(k-1))>0
                    Penalty4=Penalty4+HConstraint(1,1+i+10*(k-1));
                end
            end
        end
        
        for i=2:nNodes                                                     % 最大爬升俯冲角约束
            a_i=[x(i)-x(i-1),y(i)-y(i-1)];
            delta_z=abs(z(i)-z(i-1));
            if delta_z/norm(a_i)>1
            Penalty5=Penalty5+delta_z/norm(a_i);
            end
        end

        PopObj(1,1)=PopObj(1,1)+10^8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(1,2)=PopObj(1,2)+10^8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
        PopObj(1,3)=PopObj(1,3)+10^8*(Penalty1+Penalty2+Penalty3+Penalty4+Penalty5+Penalty6);
       
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


