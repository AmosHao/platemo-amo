classdef MOEADSS < ALGORITHM
% <multi/many> <real/binary/permutation>
% <algorithm> <M>
% Multiobjective evolutionary algorithm based on decomposition


%------------------------------- Reference --------------------------------
% Q. Zhang and H. Li, MOEA/D: A multiobjective evolutionary algorithm based
% on decomposition, IEEE Transactions on Evolutionary Computation, 2007,
% 11(6): 712-731.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATL  AB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
   
%% Parameter setting

    %% Generate the weight vectors
    [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
    T = ceil(Problem.N/10);
    W = 1./W./repmat(sum(1./W,2),1,size(W,2));
    %% Detect the neighbours of each solution
     B = pdist2(W,W);
     [~,B] = sort(B,2);
     B = B(:,1:T);
     balt=0.9;
%
    %% Generate random population
    Population = Problem.Initialization();
    Z = min(Population.objs,[],1);
    flag=0;
    flag2=0;
    P_alf=zeros(1,Problem.N);
    P_alf_flag=zeros(1,Problem.N);
    S_alf=zeros(1,Problem.N)-1;
    num_alf=zeros(1,Problem.N);
    alf=0;
    his_num=zeros(1,Problem.maxFE);
    j=0;
    mod1=0;
    mean_alf=0;
    %% Optimization
    while Algorithm.NotTerminated(Population)
     
        % For each solution
        for i = 1 : Problem.N
             alf=P_alf(i);      
             B1=B(1:Problem.N,1:size(B,2)-alf);%得到调整后的邻域T
             T1=T-alf;
              if rand < balt
                    P = B1(i,randperm(size(B1,2)));
                       %P = B1(i,1:T1);
              else
                    qwq=randperm(Problem.N);
                    P = B1(qwq(1),randperm(size(B1,2)));               
              end
               if(Problem.FE>=Problem.maxFE*0.2)
                       Offspring = OperatorGAhalf(Problem,Population(P(1:2)));
                     %Offspring = OperatorDE(Problem,Population(i),Population(P(1)),Population(P(2))); 
               else
                    Offspring = OperatorDE(Problem,Population(i),Population(P(1)),Population(P(2)));%DE算子
                    %Offspring = GAhalf(Population(P(1:2)));  
               end
          Z = min(Z,Offspring.obj);
          Zmax  = max(Population.objs,[],1);
          [~,i2]=min(max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),Problem.N,1).*W(1:Problem.N,:),[],2));
          P = B1(i2,1:T1);
          g_old = max(abs(Population(P).objs-repmat(Z,T1,1))./repmat(Zmax-Z,T1,1).*W(P,:),[],2);
          g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T1,1).*W(P,:),[],2);
                    
           %z = min([z;Offspring.objs],[],1);
           a=size(find(g_old>g_new),1);
           rate=(a./(1+exp(-10 *(Problem.FE/Problem.maxFE-0.5)) ));      
%            if(a>0&&Global.gen>=Global.maxgen*0.2)
             if(a>0&&Problem.FE>=Problem.maxFE*0.2)
               rate=(1-T1/T)*rate;
               a;
             T2=floor(rate);  
             if(T2==0)
                 T2=1;
             end
             ret=g_old-g_new;
             [~,ret2]=sort(ret);
             ret2=ret2(size(ret2)+1-size(find(ret>0),1):end);
             Population(P(ret2(1:T2))) = Offspring;

            % Population(P(find(g_old>=g_new,T2))) = Offspring;
           end
           
           
           for rep=1:T1
                if(g_old(rep)>g_new(rep))
                     Offspring_re=Population(P(rep));
                     Population(P(rep))=Offspring;
                    flag=1;
                    break;
                end        
            end  
           while(flag==1)
                Offspring=Offspring_re;
                Zmax  = max(Population.objs,[],1);
                [~,i2]=min(max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),Problem.N,1).*W(1:Problem.N,:),[],2));
                B1=B(1:Problem.N,1:size(B,2)-P_alf(i2));%得到调整后的邻域T
                T1=T-P_alf(i2);
                P = B1(i2,1:T1);
                g_old = max(abs(Population(P).objs-repmat(Z,T1,1))./repmat(Zmax-Z,T1,1).*W(P,:),[],2);
                g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T1,1).*W(P,:),[],2); 
                for rep=1:T1
                    if(g_old(rep)>g_new(rep))
                        Offspring_re=Population(P(rep));
                        Population(P(rep)) = Offspring;
                        flag=1;
                        rep=100;
                        break;
                    end        
                end
                if(rep~=100) 
                    flag=0;
                end 
           end
          if(a>0)
              P_alf_flag(i)=1;
              S_alf(i)=alf;
              num_alf(i)=a;%i有问题，与K一样
          else
               P_alf_flag(i)=0;
          end 
          if(a>0&&flag2==1)
%               S_alf(i)=alf;
%               num_alf(i)=a;%i有问题，与K一样
              if(max(num_alf)==0) 
                  max_num=1;
               else
                  max_num=max(num_alf);              
              end  
               mean_alf;
               alf2=ceil(alf+(mean_alf-alf)*((1-a/T1)));
               %alf2=ceil(alf+0.8*(mean_alf-alf)*(T1/T));
               %alf2=(alf+mean_alf)/2;
               alf=ceil(normrnd(alf2,0.2));
               while alf<0
                  alf=-alf;%向上取整
               end
               if alf>T-3
                  alf=1;
               end   
          end
          if(a==0&&flag2==1)  
               if(max(num_alf)==0) 
                  max_num=1;
               else
                  max_num=max(num_alf);              
               end  
               
               alf2=ceil(alf+(mean_alf-alf)*((1-a/T1)*a/max_num));
               %alf2=ceil(alf+0.8*(mean_alf-alf)*(T1/T));
               %alf2=(alf+mean_alf)/2;
               alf=ceil(normrnd(alf2,0.2));
               while alf<0
                  alf=-alf;
               end
               if alf>T-3
                  alf=1;
               end          
          end                     
               P_alf(i)=alf;
      end
 %---------------------------------------该种群进入下一代-----------------------------------    
        P_alf;
        j=j+1;
       % his_num(j)=sum(num_alf);
       his_num(j)=sum(P_alf_flag);
        if(j==1)
            j1=1;
        else
            j1=j-1;
        end
    %    if(sum(num_alf)<=((his_num(j1))./(1+exp(10*(Global.gen/Global.maxgen-0.5)))))
    %sum(P_alf_flag);
    
        if(mod1>0)
            mod1=mod1-1;
            flag2=0;
        end
    mod1;
      %  if(sum(P_alf_flag)<=((min(his_num(1:j)))./(1+exp(10*(Global.gen/Global.maxgen-0.5)))))
        %if(sum(P_alf_flag)<=min(his_num(1:j)))
        if(sum(P_alf_flag)<=his_num(j1)&&Problem.FE>=Problem.maxFE*0.2&&mod1==0)
            flag2=1;
            mod1=5;
        else
            flag2=0;
        end 
     
           if(sum(num_alf)>0)
           %  flag2=0;
              mean_alf=0;
%              S_alf(find(S_alf>=0));
%              num_alf(find(S_alf>=0));
              sum_val=sum(num_alf(find(S_alf>=0))./(T-S_alf(find(S_alf>=0))));
             for w=1:size(S_alf,2)
                if(S_alf(w)~=-1) 
                  % mean_alf=mean_alf+(num_alf(w)/(T- S_alf(w))/sum(num_alf/(T- S_alf))*S_alf(w));
                  mean_alf=mean_alf+(num_alf(w)/(T- S_alf(w))/sum_val*S_alf(w));
                end
             end
           end
           mean_alf;
        de=sum(P_alf_flag)/Problem.N;
        ga=1-de;
        S_alf=zeros(1,Problem.N)-1;
        num_alf=zeros(1,Problem.N);
    end
        end
    end
end

















