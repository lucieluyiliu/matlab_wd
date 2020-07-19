%this function takes jackknife output (coefficient, tstats and confidence interval)
%as input and outputs a latex table of selected group of countries.


function[]=jacknifetable(io,info,group)

  varlist={'SIZE','BM','fht','ivol','ret','roe','cash','invop','capex','ppe','rdratio','DIV','dy','close','lev','analyst','fxsale','ADR'};
if group=='dm'
    iso=info.dmiso;
    idx=info.dmidx;
else
    iso=info.emiso;
    idx=info.emidx;
end
nvar=size(io.beta,1);
ncty=size(io.beta,2);
N=io.N;
betadisp=printbeta(io);
tstatdisp=printtstat(io);
cidisp=printinterval(io);
    n=size(io.beta,1);  %number of variables   
    
 fileout = fopen(sprintf('jacknife%s %s.tex',group,io.ynames{:}),'w');
 
fprintf(fileout, '\\begin{sidewaystable}[h!] \n');
fprintf(fileout,'\\centering\n ');
fprintf(fileout,'\\resizebox{1\\textwidth}{!}{\n ');
fprintf(fileout,'\\renewcommand{\\arraystretch}{1}\n ');
fprintf(fileout,'\\begin{tabular}{l*{%d}{c}}\n',length(idx));
fprintf(fileout,'\\hline\\hline \n');
fprintf(fileout,strcat(repmat('& %s ',1,length(iso)),'\\\\ \n'),iso{:});
for i=1:nvar
fprintf(fileout,strcat('%s',repmat( ' %s ',1,length(iso)),'\\\\ \n'),varlist{i},betadisp(i,idx));
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
