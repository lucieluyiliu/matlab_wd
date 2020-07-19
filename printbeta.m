function out=printbeta(est)
beta=est.beta;
nctry=size(beta,2);
nvar=size(beta,1);
out=strings(nvar,nctry);
nselected=sum(est.selected,2);
[~,order]=sort(nselected,'descend');
for i=1:nvar
    for j=1:nctry
      
    if est.selected(order(i),j)==1
         issig=est.pvalue(order(i),j)<0.05 ; 
         if issig%if significant
         out(i,j)=sprintf('&\\textbf{%8.3f}',beta(order(i),j));
         else
             out(i,j)=sprintf('&%8.3f',beta(order(i),j));
         end
    else
        out(i,j)='&';
         
    end
end

end
end