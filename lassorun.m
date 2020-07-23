%table=readtable('holdings_fundamental_clean.csv');
%p=findpath('paneldata');
%addpath(p);
% lasobeta=zeros(nvar,nctry,7);
% lasotstat=zeros(nvar,nctry,7);
% lasop=zeros(7,nvar,nctry);
% lassoselected=boolean(zeros(nvar,nctry,7));
% lasoN=zeros(nctry,7);


  varlist={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','analyst','fxsale','ADR'};
nvar=length(varlist);
lassobeta=zeros(nvar,nctry,7);
lassotstat=zeros(nvar,nctry,7);
lassop=zeros(nvar,nctry,7);
lasssoselected=boolean(zeros(nvar,nctry,7));
lassoN=zeros(7,nctry);
yynames={'io','io_dom','io_for_us','io_for_nus','io_passive','io_active','retail'};

for i=1:7 %y regressor index

lassoout(i).ynames=yynames{i};
lassoout(i).beta=zeros(nvar,nctry);
lassoout(i).tstat=zeros(nvar,nctry);
lassoout(i).pvalue=zeros(nvar,nctry);
lassoout(i).selected=boolean(zeros(nvar,nctry));
for ci=1:length(ctrylist)
    fprintf('%s%d \n',yynames{i},ci);
    %no variable selected for hungary ci=17
    
  ctry=cell2mat(ctrylist(ci));
  ctrytbl=nustable(strcmp(nustable.ctry,ctry),:);
  year=table2array(ctrytbl(:,'year'));
  yeardummy0=dummyvar(recode(year));
  dscd=table2array(ctrytbl(:,'dscd'));
  [~,~,numid]=unique(dscd);

    

X=table2array(ctrytbl(:,varlist));

y=table2array(ctrytbl(:,yynames{i}));


%remove nan
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
%initialize lasso parameter by rule of thumb
Dt = dummyvar(categorical(id));
yeardummy=yeardummy0(idx,:);

X=X(idx,:); %if any of the columns is nan
tmpX=[X,yeardummy(:,2:end)];
Xkeep=findNonCollinear(tmpX);
tmp=tmpX(:,Xkeep);
%ZZ=ones(size(X,1),1); %partial out

ZZ=[ones(size(X,1),1),tmp(:,length(varlist)+1:end)]; %remove dummies that are colinear
%y=y-ZZ*((ZZ'*ZZ)\(ZZ'*y)); %partial out for y
%X=X-ZZ*((ZZ'*ZZ)\(ZZ'*X)); %partial out constant and year dummy
%partial out X only but should do both

n_cluster=size(Dt,2);
p=size(X,2);
n=length(y);
gamma = .1/log(n_cluster);
lambda = 1.1*2*sqrt(n)*norminv(1-gamma/(2*p));

 blasso = feasiblePostLasso(y,X,'lambda',lambda,'MaxIter',100,'UpsTol',1e-4,'ClusterVar',id);
 lassoselect = abs(blasso(1:length(varlist))) > 1e-4;

[blasso,bpost,use,bpartlasso,bpartpost,alasso,apost,Ups]=mylasso(y,X,'lambda',lambda,'UpsTol',1e-4,'MaxIter',100,'ClusterVar',id,'partial',ZZ);
%lassoselect=use;
lassoout(i).selected(:,ci)=lassoselect;
colnames=varlist(lassoselect);
y=table2array(ctrytbl(:,yynames{i}));
X=table2array(ctrytbl(:,colnames));
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
Xkeep=findNonCollinear(X);
X=X(:,Xkeep);
time=year(idx);
plasso=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);
wald=waldsigtest(plasso);
lassoout(i).beta(lassoselect,ci)=plasso.coef(1:length(colnames));
lassoout(i).tstat(lassoselect,ci)=wald.value(1:length(colnames));
lassoout(i).pvalue(lassoselect,ci)=wald.p(1:length(colnames));
lassoout(i).N(ci)=length(y);

% tmp=zeros(nvar,1);
% tmp(lassoselect)=plasso.coef(1:length(colnames));
% lassobeta(:,ci,i)=tmp;
% tmp=zeros(nvar,1);
% tmp(lassoselect)=wald.value(1:length(colnames));
% lassotstat(:,ci,i)=tmp;
% tmp=zeros(nvar,1);
% tmp(lassoselect)=wald.p(1:length(colnames));
% lassop(:,ci,i)=tmp;
% lassoN(ci,i)=size(y,1);

end

end