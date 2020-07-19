
n=length(y);
p=size(X,2);
Dt = dummyvar(categorical(id));
n_cluster=size(Dt,2);

gamma = .1/log(n_cluster);
lambda = c*2*sqrt(n)*norminv(1-gamma/(2*p));
X=[table2array(ctrytbl(:,varlist))];
y=table2array(ctrytbl(:,ynames));

idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan

Xkeep=findNonCollinear(X);
X=X(:,Xkeep);
blasso = feasiblePostLasso(y,X,'lambda',lambda,'MaxIter',100,'UpsTol',1e-4,'ClusterVar',id);
lassoselect = abs(blasso(1:length(varlist))) > 1e-4;
use=lassoselect(1:length(varlist));
varlist(use)
Xs=X(:,use);
pols=X(:,use)\y;

