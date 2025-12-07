function [PivotSol,Distance,r,t] = F_PivotSol(FunctionValue,FrontValue,MaxFront,r,t)
% This function finds all the pivot solutions in each front

[N,M] = size(FunctionValue);


PivotSol = false(1,N);

Distance = zeros(1,N);
AR = zeros(1,N);


for i = 1 : MaxFront   
    Temp = find(FrontValue==i);
    if length(Temp) <= M     
        PivotSol(Temp) = 1;
    else
        result = [];
        FF = FunctionValue(Temp,:);
        CS = 1 : M;
        for ii = 1:M
            CS1 = circshift(CS,[1,-(ii-1)]);
            PP = FF(:,CS1);
            [~,result(:,ii)] = sortrows(PP);
        end
        
        [~, AAA] = sort(result, 1);
        AR(Temp) = sum(AAA,2)';
        
        Dist  = F_dist(FF);
        Distance(Temp) = Dist';
        
        nn = size(FF,1);
        
        ARmin    = min(AR(Temp));
        ARmean   = mean(AR(Temp));
        ARmax    = max(AR(Temp));

        rate        =  ((ARmean - ARmin)/ARmax)^(ceil(M/nn));
        % update the parameter r, and calculate R
        Fmax = max(FunctionValue(Temp,:),[],1);
        Fmin = min(FunctionValue(Temp,:),[],1);
        
        if t(i) == -1
            r(i) = 1;
        else
            r(i) = r(i)/exp((1-t(i)/rate)/M);
        end
   
        R = (Fmax-Fmin).*r(i);
        
        [null,Rank] = sort(AR(Temp),'ascend');
        Choose = zeros(1,length(Rank));
        Remain = ones(1,length(Rank));
        for j = Rank
            if Remain(j)
                for k = 1 : length(Temp)
                    if abs(FunctionValue(Temp(j),:)-FunctionValue(Temp(k),:)) <= R
                        Remain(k) = 0;
                    end
                end
                Choose(j) = 1;
            end
        end
        Choose(Rank(find(Choose(Rank)==1,1,'last'))) = 0;
        PivotSol(Temp(Choose==1)) = 1;
        % update the parameter t
        t(i) = sum(Choose)/length(Temp);
    end
end
end

