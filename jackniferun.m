%prepare data
table=readtable('holdings_fundamental_clean.csv');
p=findpath('paneldata');
addpath(p);

%jacknife procedure 
nustable=table(table.isus==0,:);
ctrylist=table2array(unique(nustable(:,'ctry')));
nctry=length(ctrylist);
  varlist={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','analyst','fxsale','ADR'};
nvar=length(varlist);
yynames={'io','io_dom','io_for_us','io_for_nus','io_passive','io_active','retail'};
for i=1:7 %y regressor index
out(i).beta=zeros(nvar,nctry);
out(i).low=zeros(nvar, nctry);
out(i).high=zeros(nvar, nctry);
out(i).tstat=zeros(nvar,nctry);
out(i).p=zeros(nvar,nctry);
out(i).selected=boolean(zeros(nvar,nctry));
out(i).N=zeros(nctry,1);
for ci=1:length(ctrylist)
  ctry=cell2mat(ctrylist(ci));
  ctrytbl=nustable(strcmp(nustable.ctry,ctry),:);
  year=table2array(ctrytbl(:,'year'));
  yeardummy=dummyvar(categorical(year));
  dscd=table2array(ctrytbl(:,'dscd'));
  [~,~,numid]=unique(dscd);
out(i).ynames=yynames{i};
  [out(i).low(:,ci),out(i).high(:,ci)]=jacknife(ctrytbl, out(i).ynames);
  out(i).selected(:,ci)=~((out(i).low(:,ci)<=0)&(0<=out(i).high(:,ci)));

%arrange variable by number of countries being selected

colnames=varlist(out(i).selected(:,ci));
out(i).ynames=yynames{i};
y=table2array(ctrytbl(:,out(i).ynames));

X=[table2array(ctrytbl(:,colnames)),yeardummy(:,2+omityear:end)];
%X=[table2array(ctrytbl(:,colnames))];  %exclude one dummy
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
Xkeep=findNonCollinear(X);
X=X(:,Xkeep);
time=year(idx);
reg3=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);
wald=waldsigtest(reg3);
out(i).beta(out(i).selected(:,ci),ci)=reg3.coef(1:length(colnames));  %omit year dummy and constant
out(i).tstat(out(i).selected(:,ci),ci)=wald.value(1:length(colnames)); %if omit then total year dummy-2?
out(i).p(out(i).selected(:,ci),ci)=wald.p(1:length(colnames));
out(i).N(ci)=reg3.N;
end


end
save('jacknkife_withdummy.mat')