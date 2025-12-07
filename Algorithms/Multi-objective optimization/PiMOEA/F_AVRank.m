function SRANK = F_AVRank(FunctionValue)

[N,M] = size(FunctionValue);

CS = 1 : M;
for i = 1:M
    CS1 = circshift(CS,[1,-(i-1)]);
    PP = FunctionValue(:,CS1);
    [~, result(:,i)] = sortrows(PP);
end

[~, AAA] = sort(result, 1);

SRANK = sum(AAA,2);
    
end