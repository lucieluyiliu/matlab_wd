function out=printinterval(est,order)
low=est.low;
high=est.high;
nctry=size(low,2);
nvar=size(low,1);
out=strings(nvar,nctry);

for i=1:nvar
    for j=1:nctry
     if est.selected(order(i),j)==1 
            out(i,j)=sprintf('&[%8.3f, %8.3f ]',[low(order(i),j),high(order(i),j)]);
        else
            out(i,j)='&';
        
     end
    end
end


end