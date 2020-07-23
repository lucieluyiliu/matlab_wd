function b = feasiblePostLasso(y,x,varargin)
seed=(20200716);
[n,p] = size(x);

[Verbose,MaxIter,UpsTol,TolZero,lambda,beta0,clusterVar] = process_options(varargin,'Verbose',0,...
    'MaxIter',15,'UpsTol',1e-4,'TolZero',1e-4,'lambda',LassoSimulateLambda( x, 1, 10000, n, .05 ),...
    'beta0',[],'clusterVar',[]);

if isempty(clusterVar)
    if isempty(beta0)
        Syx = x.*(y*ones(1,p));
        Ups0 = sqrt(sum(Syx.^2)/n);
        b = LassoShooting2(x, y , .5*lambda, Ups0, 'Verbose', Verbose, 'beta', beta0);
        use = abs(b) > TolZero;
        
        e = y-x(:,use)*(x(:,use)\y);
    else
        e = y - x*beta0;
        Syx = x.*(e*ones(1,p));
        Ups0 = sqrt(sum(Syx.^2)/n);

        b = LassoShooting2(x, y , lambda, Ups0, 'Verbose', Verbose, 'beta', beta0);
        %use = abs(b) > TolZero;
        use = abs(b) > TolZero;
        e = y-x(:,use)*(x(:,use)\y);
    end
    
    kk = 1;
    Syx = x.*(e*ones(1,p));
    Ups1 = sqrt(sum(Syx.^2)/n);
    while norm(Ups0-Ups1) > UpsTol && kk < MaxIter
       
        d0 = norm(Ups0-Ups1);
        b = LassoShooting2(x, y , lambda, Ups1 , 'Verbose', Verbose, 'beta', b);
        use = abs(b) > TolZero;
        
        e = y-x(:,use)*(x(:,use)\y);
        Ups0 = Ups1;
        Syx = x.*(e*ones(1,p));
        Ups1 = sqrt(sum(Syx.^2)/n);
        kk = kk+1;
     
        d1 = norm(Ups0 - Ups1);
        if d1 == d0
            Ups0 = Ups1;
        end
        
    end    
else
    Dt = dummyvar(categorical(clusterVar));
    n_clust=size(Dt,2);
    gamma=0.1/log(n_clust);
    p=size(x,2);
     n=size(x,1);
    lambda=2*1.1*sqrt(n)*norminv(1-gamma/(2*p));
    if isempty(beta0)
        Syx = x.*(y*ones(1,p));
        Ups0 = sqrt(sum((Dt'*Syx).^2)/n);
        b = LassoShooting2(x, y , .5*lambda, Ups0, 'Verbose', Verbose, 'beta', beta0);
        use = abs(b) > TolZero;
        e = y-x(:,use)*(x(:,use)\y);
    else
        e = y - x*beta0;
        Syx = x.*(e*ones(1,p));
        Ups0 = sqrt(sum((Dt'*Syx).^2)/n);
        b = LassoShooting2(x, y , lambda, Ups0, 'Verbose', Verbose, 'beta', beta0);
        use = abs(b) > TolZero;
        e = y-x(:,use)*(x(:,use)\y);
    end
    
    kk = 1;
    Syx = x.*(e*ones(1,p));
    Ups1 = sqrt(sum((Dt'*Syx).^2)/n);
    while norm(Ups0-Ups1) > UpsTol && kk < MaxIter
        
        d0 = norm(Ups0-Ups1);
        b = LassoShooting2(x, y , lambda, Ups1 , 'Verbose', Verbose, 'beta', b);
        use = abs(b) > TolZero;
        e = y-x(:,use)*(x(:,use)\y);
        Ups0 = Ups1;
        Syx = x.*(e*ones(1,p));
        Ups1 = sqrt(sum((Dt'*Syx).^2)/n);

        kk = kk+1;
        d1 = norm(Ups0 - Ups1);
        if d1 == d0
            Ups0 = Ups1;
        end
        
    end
    
end
end