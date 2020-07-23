%this script take lasso and variable selection output as given and produce
%latex tables

for i=1:7

 jacknifetable(jkout(i),'dm');
 jacknifetable(jkout(i),'em');   
end

for i=1:7
 lassotable(lassoout(i),'dm');
 lassotable(lassoout(i),'em');
end



for i=1:7
 selectiontable(jkout(i),lassoout(i),'dm');
 selectiontable(jkout(i),lassoout(i),'em');
end

 

