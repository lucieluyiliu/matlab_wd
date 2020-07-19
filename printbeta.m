function out=printbeta(io)
beta=io.beta;
nctry=size(beta,2);
nvar=size(beta,1);
out=strings(nvar,nctry);
for i=1:nvar
    for j=1:nctry
      
    if io.selected(i,j)==1
         issig=io.p(i,j)<0.05; 
         if issig%if significant
         out(i,j)=sprintf('&\\textbf{%8.3f}',beta(i,j));
         else
             out(i,j)=sprintf('&%8.3f',beta(i,j));
         end
    else
        out(i,j)='&'
         
    end
end

end
end