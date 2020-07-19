
%prepare data
table=readtable('holdings_fundamental_clean.csv');


%focus on non-US firms
%*****************************
%prepare data

%*****************************
nustable=table(table.isus==0,:);

dscd=table2array(nustable(:,'dscd'));
[~,~,numid]=unique(dscd);
io=table2array(nustable(:,'io'));
io_dom=table2array(nustable(:,'io_dom'));
io_for_us=table2array(nustable(:,'io_for_us'));
io_for_nus=table2array(nustable(:,'io_for_nus'));
io_active=table2array(nustable(:,'io_active'));
io_passive=table2array(nustable(:,'io_passive'));
retail=table2array(nustable(:,'retail'));





%prepare year and country fixed effects

year=table2array(nustable(:,'year'));
country=table2array(nustable(:,'ctry'));

nustable(:,'ctry')=strrep(table2array(nustable(:,'ctry')),' ','');
yeardummy=dummyvar(categorical(year));
yearlist=unique(year);
ctrylist=unique(country);
yeardummylist=cell(1,length(yearlist));

for i=1:length(yearlist)
    yeardummylist{i}=sprintf('year=%d',yearlist(i));
end
countrydummy=dummyvar(categorical(country));
ctrydummylist=cell(1,length(ctrylist));
for i=1:length(ctrylist)
    ctrydummylist{i}=sprintf('%s',ctrylist{i});
end

%all countries io regression
%1. country level variable
ynames={'io'};
y=table2array(nustable(:,ynames));
colnames={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','ADR','analyst','fxsale','legal','disc','distance','english','gdp','mktcap','investable','inflation','pol'};
X=[table2array(nustable(:,colnames))];  %exclude one dummy
xnames=[colnames,yeardummylist{2:end}];
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
time=year(idx);
io1=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);

ynames={'io'};

io1.yname=ynames;
io1.xnames=xnames;
estdisp(io1)



%2. country fixed effects

y=io;
colnames={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','ADR','analyst','fxsale'};
X=[table2array(nustable(:,colnames)),yeardummy(:,2:end),countrydummy(:,2:end)];  %exclude one dummy
ynames={'retail'};
xnames=[colnames,yeardummylist{2:end},ctrydummylist(2:end)];
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
time=year(idx);
io2=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);



io2.yname=ynames;
io2.xnames=xnames;
estdisp(io2)


%jacknife procedure 

ctrylist=table2array(unique(nustable(:,'ctry')));
nctry=length(ctrylist);
  varlist={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','analyst','fxsale','ADR'};
nvar=length(varlist);
beta=zeros(nvar,nctry);
low=zeros(nvar, nctry);
high=zeros(nvar, nctry);
tstat=zeros(nvar,nctry);
p=zeros(nvar,nctry);
selected=boolean(zeros(nvar,nctry));
N=zeros(nctry,1);


for ci=1:length(ctrylist)
  ctry=cell2mat(ctrylist(ci));
%   if strcmp(ctry,'COLOMBIA' )
%     omityear=1;
% else omityear=0;
%end
ynames='io';
  ctrytbl=nustable(strcmp(nustable.ctry,ctry),:);
  year=table2array(ctrytbl(:,'year'));
  yeardummy=dummyvar(recode(year));
  dscd=table2array(ctrytbl(:,'dscd'));
  [~,~,numid]=unique(dscd);

  [low(:,ci),high(:,ci)]=jacknife(ctrytbl,ynames);
  selected(:,ci)=~((low(:,ci)<=0)&(0<=high(:,ci)));

%arrange variable by number of countries being selected

colnames=varlist(selected(:,ci));
ynames={'io'};
y=table2array(ctrytbl(:,ynames));
X=[table2array(ctrytbl(:,varlist)),yeardummy(:,2:end)];
idx=~any(isnan(X), 2); %valid non-nan index
y=y(idx);
id=numid(idx);
X=X(idx,:); %if any of the columns is nan
Xkeep=findNonCollinear(X);
X=X(:,Xkeep);
time=year(idx);
reg3=panel(id,time,y,X,'po','vartype','cluster','clusterid',id);
wald=waldsigtest(reg3);
beta(selected(:,ci),ci)=reg3.coef(1:length(colnames));  %omit year dummy and constant
tstat(selected(:,ci),ci)=wald.value(1:length(colnames)); %if omit then total year dummy-2?
p(selected(:,ci),ci)=wald.p(1:length(colnames));
N(ci)=reg3.N;
end

info.developed=["AUSTRALIA","AUSTRIA", "BELGIUM","CANADA","DENMARK","FINLAND","FRANCE","GERMANY","HONGKONG","IRELAND","ISRAEL","ITALY","JAPAN","LUXEMBOURG","NETHERLANDS,""NEWZEALAND","NORWAY","PORTUGAL","SINGAPORE","SPAIN","SWEDEN","SWITZERLAND","UNITEDKINGDOM"];
info.emerging={'BRAZIL','CHILE','CHINA','COLOMBIA','CZECHREPUBLIC','EGYPT','GREECE','HUNGARY','INDIA','INDONESIA','MALAYSIA','MEXICO','MOROCCO','PERU','PHILIPPINES','POLAND','RUSSIA','SOUTHAFRICA','SOUTHKOREA','TAIWAN','THAILAND','TURKEY'};
info.developediso={'AU','AT','BE','CA','DK','FI','FR','DE','HK','IE','IL','IT','JP','LU','NL','NZ','NO','PT','SG','ES','SE','CH','UK'};
info.emergingiso={'BR','CL','CN','CO','CZ','EG','GR','HU','IN','ID','MY','MX','MA','PE','PH','PL','RU','ZA','KR','TW','TH','TR'};
info.dmidx=zeros(length(developed),1);
info.emidx=zeros(length(emerging),1);
for i=1:length(developed)
dmidx(i)=find(strcmp(ctrylist,developed(i)));
end
for i=1:length(emerging)
emidx(i)=find(strcmp(ctrylist,emerging(i)));
end

