%this function takes jackknife output (coefficient, tstats and confidence interval)
%as input and outputs a latex table of selected group of countries.


function[]=jacknifetable(est,group)
ctrylist={'AUSTRALIA','AUSTRIA','BELGIUM' ,'BRAZIL' ...
    'CANADA'  ,'CHILE'  ,'CHINA','COLOMBIA','CZECHREPUBLIC',...
   'DENMARK','EGYPT','FINLAND','FRANCE','GERMANY','GREECE' ,...
   'HONGKONG','HUNGARY','INDIA','INDONESIA', 'IRELAND'  ,...
    'ISRAEL' , 'ITALY', 'JAPAN', 'LUXEMBOURG', 'MALAYSIA' ...
    'MEXICO','MOROCCO','NETHERLANDS' , 'NEWZEALAND', 'NORWAY' ,...
    'PERU', 'PHILIPPINES', 'POLAND', 'PORTUGAL', 'RUSSIA' ,...
    'SINGAPORE' , 'SOUTHAFRICA' , 'SOUTHKOREA', 'SPAIN',...
    'SWEDEN' , 'SWITZERLAND', 'TAIWAN',  'THAILAND', 'TURKEY',...
    'UNITEDKINGDOM'};
developed={'AUSTRALIA','AUSTRIA','BELGIUM',...
'CANADA','DENMARK','FINLAND','FRANCE','GERMANY','HONGKONG','IRELAND','ISRAEL','ITALY','JAPAN','LUXEMBOURG','NETHERLANDS','NEWZEALAND','NORWAY','PORTUGAL','SINGAPORE','SPAIN','SWEDEN','SWITZERLAND','UNITEDKINGDOM'};
emerging={'BRAZIL','CHILE','CHINA','COLOMBIA','CZECHREPUBLIC','EGYPT','GREECE','HUNGARY','INDIA','INDONESIA','MALAYSIA','MEXICO','MOROCCO','PERU','PHILIPPINES','POLAND','RUSSIA','SOUTHAFRICA','SOUTHKOREA','TAIWAN','THAILAND','TURKEY'};
dmiso={'AU','AT','BE','CA','DK','FI','FR','DE','HK','IE','IL','IT','JP','LU','NL','NZ','NO','PT','SG','ES','SE','CH','UK'};
emiso={'BR','CL','CN','CO','CZ','EG','GR','HU','IN','ID','MY','MX','MA','PE','PH','PL','RU','ZA','KR','TW','TH','TR'};
dmname={'Australia','Austria','Belgium','Canada','Denmark','Finland','France','Germany','HK','Ireland','Israel','Italy','Japan','Luxembourg','Netherlands','New Zealand','Norway','Portugal','Singapore','Spain','Sweden','Switzerland','UK'};
emname={'Brazil','Chile','China','Clombia','Czech Republic','Egypt','Greece','Hungary','India','Indonesia','Malaysia','Mexico','Morocco','Peru','Philippines','Poland','Ruassia','South Africa','South Korea','Taiwan','Thailand','Turkey'};

dmidx=zeros(length(developed),1);
emidx=zeros(length(emerging),1);
 %sort variable by number of markets selected
for i=1:length(developed)
dmidx(i)=find(strcmp(ctrylist,developed(i)));
end
for i=1:length(emerging)
emidx(i)=find(strcmp(ctrylist,emerging(i)));
end
  varlist={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','analyst','fxsale','ADR'};
if group=='dm'
    iso=dmiso;
    idx=dmidx;
    name=dmname;
else
    iso=emiso;
    idx=emidx;
    name=emname;
end
nvar=size(est.beta,1);
ncty=size(est.beta,2);
N=est.N;
  nselected=sum(est.selected(:,idx),2);
[~,order]=sort(nselected,'descend');  
betadisp=printbeta(est,order);
tstatdisp=printtstat(est,order);
cidisp=printinterval(est,order);
    n=size(est.beta,1);  %number of variables   

 fileout = fopen(sprintf('jacknife%s %s.tex',group,est.ynames),'w');
 
fprintf(fileout, '\\begin{sidewaystable}[h!] \n');
fprintf(fileout,'\\centering\n ');
fprintf(fileout,'\\resizebox{1\\textwidth}{!}{\n ');
fprintf(fileout,'\\renewcommand{\\arraystretch}{1}\n ');
fprintf(fileout,'\\begin{tabular}{l*{%d}{c}}\n',length(idx));
fprintf(fileout,'\\hline\\hline \n');
fprintf(fileout,strcat(repmat('& %s ',1,length(iso)),'\\\\ \n'),name{:});
for i=1:nvar
fprintf(fileout,strcat('%s',repmat( ' %s ',1,length(iso)),'\\\\ \n'),varlist{order(i)},betadisp(i,idx));
fprintf(fileout,strcat(repmat('%s ',1,length(iso)),'\\\\ \n'),tstatdisp(i,idx));
fprintf(fileout, strcat(repmat('%s ',1,length(iso)),'\\\\ \n'),cidisp(i,idx));
end
fprintf(fileout, '\\hline \n');
fprintf(fileout, strcat('N',repmat('& %d ',1,length(iso)),'\\\\ \n'),N(idx));
fprintf(fileout, '\\hline\\hline \n');
%fprintf(strcat(repmat('& %s%8.3f%s ',1,length(idx)),'\\\\ \n'),[leftbf{:};io.beta(1,idx);rightbf{:}])  %beta estimation
%fprintf(strcat(repmat('& (%8.3f) ',1,length(idx)),'\\\\ \n'),io.tstat(1,idx))  %tstat 
%fprintf(strcat(repmat('& [%8.3f,%8.3f]',1,length(idx)),'\\\\ \n'),[io.low(1,idx);io.high(1,idx)]);
fprintf(fileout, '\\end{tabular}}\n');
fprintf(fileout, '\\end{sidewaystable}');
fclose(fileout);

%leftbf=arrayfun(@(x)textboldl(x),issig,'UniformOutput',false)
%fightbf=arrayfun(@(x)textboldr,issig,'UniformOutput',false)
