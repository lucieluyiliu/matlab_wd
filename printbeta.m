function out=printbeta(est,order)
betacoef=est.beta;
nctry=size(betacoef,2);
nvar=size(betacoef,1);
out=strings(nvar,nctry);
%print all country table, use order in jacknife table(dm or em)
for i=1:nvar
    for j=1:nctry
      
    if est.selected(order(i),j)==1
         issig=est.pvalue(order(i),j)<0.05 ; 
         if issig%if significant
         out(i,j)=sprintf('&\\textbf{%8.3f}',betacoef(order(i),j));
         else
             out(i,j)=sprintf('&%8.3f',betacoef(order(i),j));
         end
    else
        out(i,j)='&';
         
    end
end

end
end