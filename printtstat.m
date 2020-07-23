function out=printtstat(est,order)
tstat=est.tstat;
nctry=size(tstat,2);
nvar=size(tstat,1);
out=strings(nvar,nctry);

for i=1:nvar
    for j=1:nctry
        
        if est.selected(order(i),j)==1
            out(i,j)=sprintf('&(%8.3f)',tstat(order(i),j));
        else
            out(i,j)='&';
        end
    end
end
end
