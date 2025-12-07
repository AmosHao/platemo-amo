function testPopulationVariation(objDim,gen,pop,values,algName)
% Plot the variation figures of PF and PS
if objDim==2
    subplot(1,2,1);
    plot(values(:,1),values(:,2),'ro');
    str=sprintf('gen= %d',gen);
    title(str);
    legend(algName)
    box on;
    drawnow;
    subplot(1,2,2);
    plot3(pop(:,1),pop(:,2),pop(:,3),'bo');
    legend('obtainedPS')
    str = sprintf('gen= %d',gen);
    title(str);
    box on;
    drawnow;
elseif objDim == 3
    subplot(1,2,1);
    plot3(values(:,1),values(:,2),values(:,3),'ro');
    legend(algName)
    str = sprintf('gen= %d',gen);
    title(str);
    box on;
    drawnow;
    subplot(1,2,2);
    plot3(pop(:,1),pop(:,2),pop(:,3),'bo');
    legend('obtainedPS')
    str = sprintf('gen= %d',gen);
    title(str);
    box on;
    drawnow
end
end