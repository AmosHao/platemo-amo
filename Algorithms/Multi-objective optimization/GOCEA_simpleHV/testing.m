clus = zeros(100,2);
PopulationDec = zeros(100,12);
clus(:,1)=(1:100);
for i =1:100
    clus(i,2) = num2cell(PopulationDec(i,:));
end

disp(clus)