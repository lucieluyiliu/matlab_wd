
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
X=[table2array(nustable(:,colnames)),yeardummy(:,2:end)];  %exclude one dummy
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
ynames={'io'};
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
ctrylost=unique(nustable(:,'ctry'));
for i=1:length(ctrylist)
    

end
