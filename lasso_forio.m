
%prepare data
table=readtable('holdings_fundamental_clean.csv');
nustable(:,'ctry')=strrep(table2array(nustable(:,'ctry')),' ','');
country=table2array(nustable(:,'ctry'));
ctrylist=unique(country);
nvar=length(varlist);
beta=zeros(nvar,nctry);
low=zeros(nvar, nctry);
high=zeros(nvar, nctry);
tstat=zeros(nvar,nctry);
pvalue=zeros(nvar,nctry);
selected=boolean(zeros(nvar,nctry));
N=zeros(nctry,1);
yynames={'io','io_dom','io_for_us','io_for_nus','io_passive','io_active','retail'};
  varlist={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','analyst','fxsale','ADR'};

  
 for ci=1:length(ctrylist)
     ctry=cell2mat(ctrylist(ci));
ynames={'io'};
  ctrytbl=nustable(strcmp(nustable.ctry,ctry),:);
  year=table2array(ctrytbl(:,'year'));
  yeardummy=dummyvar(recode(year));
  dscd=table2array(ctrytbl(:,'dscd'));
  [~,~,numid]=unique(dscd);
n=length(y);
p=size(X,2);

X=table2array(ctrytbl(:,varlist));
y=table2array(ctrytbl(:,ynames));

%remove nan
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
%initialize lasso parameter by rule of thumb
Dt = dummyvar(categorical(id));
yeardummy=yeardummy(idx,:);
X=X(idx,:); %if any of the columns is nan
const=[ones(size(X,1),1),yeardummy(:,2:end)]; %partial out
y=y-mean(y);
X=X-const*((const'*const)\(const'*X)); %partial out constant and year dummy
%remove colinear
n_cluster=size(Dt,2);
gamma = .1/log(n_cluster);
lambda = 1.1*2*sqrt(n)*norminv(1-gamma/(2*p));

blasso = feasiblePostLasso(y,X,'lambda',lambda,'MaxIter',100,'UpsTol',1e-4,'ClusterVar',id);
lassoselect = abs(blasso(1:length(varlist))) > 1e-4;

selected(:,ci)=lassoselect;
colnames=varlist(selected(:,ci));
y=table2array(ctrytbl(:,ynames));
X=[table2array(ctrytbl(:,colnames)),yeardummy(:,2:end)];
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
Xkeep=findNonCollinear(X);
X=X(:,Xkeep);
time=year(idx);
plasso=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);
wald=waldsigtest(plasso);
beta(selected(:,ci),ci)=plasso.coef(1:length(colnames));
tstat(selected(:,ci),ci)=wald.value(1:length(colnames));
pvalue(selected(:,ci),ci)=wald.p(1:length(colnames));

N(ci)=plasso.N;
end