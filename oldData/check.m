%读取cubes数据
             % 文件名
            filename = 'city_environment_data_new14.xlsx';
            % 读取长方体数据
            buildings = xlsread(filename, 'Buildings');
            numBuildings = size(buildings, 1);
            % 读取圆柱体数据
            cylinders = xlsread(filename, 'Cylinders');
            numCylinders = size(cylinders, 1);
            % 读取四棱锥数据
            pyramids = xlsread(filename, 'Pyramids');
            numPyramids = size(pyramids, 1);
            % 读取球体数据
            spheres = xlsread(filename, 'Spheres');
            numSpheres = size(spheres, 1);
            % 读取禁飞区数据
            jinfei = xlsread(filename, 'jinfei');
            numjinfei = size(jinfei, 1);

dot1=[500,500,800];
dot2=[2800,4050,800];
filename = '1102_dis1o1_4_maxhjd_dot1_1_dot2_6.xlsx';
obj = xlsread(filename, 'obj_data');
dec = xlsread(filename, 'dec_data');
varDim=size(dec,2);
validtrait = xlsread(filename, 'validtrait_data');
obj_valid=obj(find(validtrait(:,1)==1),:);
dec_valid=dec(find(validtrait(:,1)==1),:);
[minValue, minIndex] = min(obj_valid(:, 3));%obj第一列最小值
[validPoints]=checknode_singlepath_close_new(dec_valid(minIndex,:),varDim,buildings,numBuildings,cylinders,numCylinders,spheres,numSpheres,pyramids,numPyramids,jinfei,numjinfei,dot1,dot2);%对当前解包含的航迹点进行碰撞检测
    a=[];b=[];              
if any(~validPoints) 
                   b=[b;0];
                   % Population(i).validtrait = 0; 
                   else   
                   b=[b;1];
                   % Population(i).validtrait = 1;
                   end