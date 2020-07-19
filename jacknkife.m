function out=jacknife(data,candidate)

%this function runs the jacknife procedure
seed=20200709;
rng(seed);
M=1000;
betas=zeros(M,1);
ctry='AUSTRIA';
ctrytbl=nustable(strcmp(nustable.ctry,ctry),:);
varlist={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','ADR','analyst','fxsale'};

dscd=table2array(ctrytbl(:,'dscd'));
[~,~,numid]=unique(dscd);
candidate='SIZE';
position=find(strcmp(varlist,candidate));

n=length(varlist);

for i=1:M
k=floor(8+rand()*(n-8));
iscontrol=(1:n)~=position;
controlvar=varlist(iscontrol);
control=datasample(controlvar,k,'Replace',false);  %sample controls

year=table2array(ctrytbl(:,'year'));
yeardummy=dummyvar(categorical(year));
yearlist=unique(year);

ynames={'io'};
colnames=[candidate,control];

y=table2array(ctrytbl(:,ynames));
X=[table2array(ctrytbl(:,colnames)),yeardummy(:,2:end)];  %exclude one dummy
xnames=[candidate,control];
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
time=year(idx);
io1=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);
betas(i)=io1.coef(1);  %candidate variable coef.
t=io1.coef/io1.stderr;
keep=t>1;
control1=varlist(keep);
colnames1=[candidate,control1];
X=[table2array(ctrytbl(:,colnames1)),yeardummy(:,2:end)];  %exclude one dummy
xnames=[candidate,control1];
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
time=year(idx);
io2=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);
betas(i)=io2.coef(1);
end

