function xtrial=testGuassianSample(mu,sigma,bounds,varDim)
    % Sample a new solution
    xtrial  = mvnrnd(mu, sigma);
    % boundary check and repair
    xupp    = bounds(2,:);
    xlow    = bounds(1,:);
    rnds    = rand(1,varDim);
    pos     = xtrial < xlow;
    xtmp    = mu - rnds.*(mu-xlow);
    xtrial(pos) = xtmp(pos);
    pos     = xtrial > xupp;
    xtmp    = mu + rnds.*(xupp-mu);
    xtrial(pos) = xtmp(pos);
end