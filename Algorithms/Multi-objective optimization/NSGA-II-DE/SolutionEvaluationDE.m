 function objvs=SolutionEvaluationDE(prob,pop)
        % Author: Zhang Hu, Harbin Institute of Technology
        % Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

            [popSize,varDim]=size(pop);
            objDim=prob.M;
            objvs=zeros(popSize,objDim);
            for i=1:popSize
                v=prob.CalObj(pop(i,1:varDim));
                objvs(i,1:objDim)=v;
            end
        end