classdef AMEA < ALGORITHM
% <multi> <real/integer>
% GOCEA
% F ---   0.5 --- The F
% CR    --- 1 --- The CR
% hisLen ---   5 --- The F
% beta ---   0.9 --- The F
% alpha    --- 0.01 --- The CR
    methods
        function main(Algorithm, Problem)
            %% parameter setting
            [F,CR,hisLen,beta,alpha] = Algorithm.ParameterSet(0.5,1,5,0.9,0.01);
            popSize=Problem.N;
            maxGens=Problem.maxFE/Problem.N/2;
            objDim=Problem.M;
            varDim=Problem.D;
            bounds=[Problem.lower;Problem.upper];
            parents=zeros(3,varDim);                 % Define the parents for DE operation
            pm=1.0/varDim;                           % Mutation probability
            %% Generate random population
            Population = Problem.Initialization();
            pop=Population.decs;
            auxPop=Population.decs;                             % Set an auxiliary population for environmental selection
            objVals=Population.objs;
            auxVals=Population.objs;
            auxPop_b=Population.decs;                            % Set an auxiliary population for beta update
            auxVals_b=Population.objs;
            tmpMatrix=zeros(4,maxGens);              % Temporary matrix for beta update
            betaSet=zeros(maxGens+1,1);              % Track the beta update
            betaSet(1,1)=beta;
            clustTime=0;                             % Record the number of performing clustering
            diffSize=popSize;                        % Record the number of solutions that survive into the next generation
           %% Optimization
           while Algorithm.NotTerminated(Population)
               if diffSize/popSize>alpha
                    sqrtPopSize=floor(sqrt(popSize)); nClust=randsample(sqrtPopSize,1);
                    while nClust==1
                        nClust=randsample(sqrtPopSize,1);
                    end
                    clustIdx=kmeans(pop,nClust,'EmptyAction','drop');
                    clustVals=objVals;
                    tmpParm=unique(clustIdx);
                    if length(tmpParm)<2            % To avoid that just one cluster exists
                        tmpParm1=tmpParm+1;
                        if tmpParm1>nClust
                            tmpParm1=tmpParm-1;
                        end
                        clustIdx(randsample(popSize,floor(popSize/2)),1)=tmpParm1;
                    end
                    clustTime=clustTime+1;
                end
                % Calculate the covariance matrices for each cluster
                rndInds=zeros(2,varDim);
                sigmaSet=cell(nClust,1);
                for i=1:nClust
                    mark=clustIdx==i;
                    lenMark=sum(mark);
                    if lenMark>1
                        neigPop=pop(mark,1:varDim);
                        sigmaSet(i,1)={cov(neigPop)};
                        if sum(rndInds(1,1:varDim))==0
                            pos=randsample(lenMark,1); rndInds(1,1:varDim)=neigPop(pos,1:varDim);
                        end
                    end
                end
                pos=randsample(popSize,1); rndInds(2,1:varDim)=pop(pos,1:varDim);
                while ismember(rndInds(2,1:varDim),rndInds(1,1:varDim),'rows')
                    pos=randsample(popSize,1); rndInds(2,1:varDim)=pop(pos,1:varDim);  % Ensure that there are at least two solutions coming from different sources
                end
                
                flg=zeros(popSize,1);
                for i=1:popSize
                    % Determine the neighbouring solutions for the current solution
                    currentSol=pop(i,1:varDim);
                    clustLabel=clustIdx(i);
                    sigma=sigmaSet{clustLabel,1};
                    
                    if ismember(rndInds(1,1:varDim),currentSol,'rows')
                        trialSol=testGuassianSample(currentSol,sigma,bounds,varDim);
                        flg(i,1)=1;
                    elseif ismember(rndInds(2,1:varDim),currentSol,'rows')
                        idx=randsample(popSize,2);
                        parents(1:2,1:varDim)=pop(idx,1:varDim);
                        parents(3,1:varDim)=currentSol;
                        trialSol=testDECrossover(parents,bounds,F,CR);
                        flg(i,1)=2;
                    elseif ~isempty(sigma)
                        if rand<beta
                            trialSol=testGuassianSample(currentSol,sigma,bounds,varDim);
                            flg(i,1)=1;
                        else
                            idx=randsample(popSize,2);
                            parents(1:2,1:varDim)=pop(idx,1:varDim);
                            parents(3,1:varDim)=currentSol;
                            trialSol=testDECrossover(parents,bounds,F,CR);
                            flg(i,1)=2;
                        end
                    else
                        idx=randsample(popSize,2);
                        parents(1:2,1:varDim)=pop(idx,1:varDim);
                        parents(3,1:varDim)=currentSol;
                        trialSol=testDECrossover(parents,bounds,F,CR);
                        flg(i,1)=2;
                    end
                    trialSol=testPolynomialMutation(trialSol,bounds,pm);
                    trialVal_total=Problem.Evaluation(trialSol);
                    trialVal=trialVal_total.objs;
                    auxPop_b(i,1:varDim)=trialSol;
                    auxVals_b(i,1:objDim)=trialVal;
                    [auxPop,auxVals]=testOriginalSMS([auxPop;trialSol],[auxVals;trialVal],popSize+1,objDim);
                end
                
                % To find out the solutions in auxiliary population auxPop that are different
                % from those in the clustering population clustPop
                diffVals=setdiff(auxVals,clustVals,'rows');
                diffSize=size(diffVals,1);
                
                % Assign the population to the next generation
                pop=auxPop;
                objVals=auxVals;
                Population=Problem.Evaluation(pop);
                
                % Update the beta
            
                auxFit=SPEA2Fitness(auxVals_b); % Calculate the raw fitness value
                tmpMatrix(1,gen)=sum(flg==1);
                tmpMatrix(2,gen)=sum(auxFit(flg==1,1));
                tmpMatrix(3,gen)=sum(flg==2);
                tmpMatrix(4,gen)=sum(auxFit(flg==2,1));
                epsilon=10^-10;
                if gen<hisLen
                    u1=sum(tmpMatrix(2,1:gen))/sum(tmpMatrix(1,1:gen));
                    u2=sum(tmpMatrix(4,1:gen))/sum(tmpMatrix(3,1:gen));
                    beta=(u2+epsilon)/(u1+u2+epsilon);
                else
                    u1=sum(tmpMatrix(2,gen-hisLen+1:gen))/sum(tmpMatrix(1,gen-hisLen+1:gen));
                    u2=sum(tmpMatrix(4,gen-hisLen+1:gen))/sum(tmpMatrix(3,gen-hisLen+1:gen));
                    beta=(u2+epsilon)/(u1+u2+epsilon);
                end
                betaSet(gen+1,1)=beta;
                % Plot the variation figures of PF and PS------------------------------
                % PopulationVariation(objDim,gen,pop,objVals,'AMEA'); 
           end
        end
    end
end