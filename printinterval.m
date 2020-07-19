function out=printinterval(io)
low=io.low;
high=io.high;
nctry=size(low,2);
nvar=size(low,1);
out=strings(nvar,nctry);
for i=1:nvar
    for j=1:nctry
     if io.selected(i,j)==1
            out(i,j)=sprintf('&[%8.3f, %8.3f ]',[low(i,j),high(i,j)]);
        else
            out(i,j)='&';
        
     end
    end
end


end