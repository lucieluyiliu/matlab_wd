function keep = findNonCollinear(Z,tol)

%choose variables to keep that are not collinear
%here I mademodification to takr into account colinear with consttant
%so that it is the same as stata
if nargin < 2
    tol = 1e-16;
end
if isempty(tol)
    tol = 1e-16;
end
n=size(Z,1);
p = size(Z,2);
keep = (1:p)';

for ii = 1:p
    use = keep ~= ii; %from first one
    ZZ=[ones(n,1),Z(:,keep(use))];
    e = Z(:,ii)-ZZ*pinv(ZZ'*ZZ)*ZZ'*Z(:,ii);
    %e = Z(:,p-ii+1)-ZZ*pinv(ZZ'*ZZ)*ZZ'*Z(:,p-ii+1);
    if sum(e.^2) < tol
        keep = setdiff(keep,ii);
        %keep = setdiff(keep,p-ii+1);
    end
   
end
