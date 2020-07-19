function [low, high]=jacknife(ctrytbl,ynames)

%this function runs the jacknife procedure
seed=20200709;
rng(seed);
M=10;

%avoid adr being 0 in a country for all firms
  varlist={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','analyst','fxsale','ADR'}; 
if sum(ctrytbl.ADR)==0
    omitadr=1;  
else
    omitadr=0;
end   
% ctry=cell2mat(ctrytbl.ctry(1));
% if strcmp(ctry,'COLOMBIA' )
%     omityear=1;
% else omityear=0;
% end
n=length(varlist)-omitadr;
low=zeros(length(varlist),1);
high=zeros(length(varlist),1);
parfor j=1:n
betas=zeros(M,1);
dscd=table2array(ctrytbl(:,'dscd'));
[~,~,numid]=unique(dscd);

candidate=varlist(j);
varpool=varlist(1:end-omitadr);  %work with varpool rather than varlist
position=find(strcmp(varpool,candidate));


for i=1:M
k=floor(8+rand()*(n-8));  

iscontrol=((1:n)~=position);
 %if omit adr then last one in varlist is omitted
controlvar=varpool(iscontrol);
control=datasample(controlvar,k,'Replace',false);  %sample controls

year=table2array(ctrytbl(:,'year'));
yeardummy=dummyvar(recode(year));

colnames=[candidate,control];


y=table2array(ctrytbl(:,ynames));
%X=[table2array(ctrytbl(:,colnames))];
%exclude one dummy
X=[table2array(ctrytbl(:,colnames)),yeardummy(:,2:end)]; 
xnames=[candidate,control];
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
Xkeep=findNonCollinear(X);
X=X(:,Xkeep);
time=year(idx);
reg1=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);
betas(i)=reg1.coef(1);  %candidate variable coef.
wald=waldsigtest(reg1);
tstat=wald.value;
keep=abs(tstat)>1;
control1=control(keep(2:length(control))); % omitcandidate and  year dummy and constant
colnames1=[candidate,control1];
X=[table2array(ctrytbl(:,colnames1)),yeardummy(:,2:end)];  %exclude one dummy

idx=~any(isnan(X), 2); %valid non-nan index
y=table2array(ctrytbl(:,ynames));
%X=[table2array(ctrytbl(:,colnames1))];
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
time=year(idx);
Xkeep=findNonCollinear(X);
X=X(:,Xkeep);
reg2=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);
betas(i)=reg2.coef(1);
end

low(j)=prctile(betas,5);
high(j)=prctile(betas,95);

   

end
if omitadr==1
     low(end)=0;
     high(end)=0;
end
end

