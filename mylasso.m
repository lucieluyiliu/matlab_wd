function [blasso,bpost,use,bpartlasso,bpartpost,alasso,apost,Ups] = mylasso(y,x,varargin)
% [blasso,bpost,use,bpartlasso,bpartpost,alasso,apost,Ups] = mylasso(y,x,varargin)
% Returns lasso (blasso) and post-lasso (bpost) estimates of coefficients
% in a linear model.  bpartlasso and bpartpost are lasso and post-lasso 
% coefficients on any variables that are partialed out before applying 
% lasso.  alasso and apost are the lasso and post-lasso intercepts when 
% 'demean' or 'standard' are used and empty otherwise.  use is a logical 
% vector for active set of x.  Ups are the estimated penalty loadings.
%
% Options (default):
% 'Verbose' (0) - controls output display during lasso updates.  0 being no
%   output and 1 and 2 corresponding to more output displayed
%
% 'MaxIter' (1) - number of lasso iterations to update penalty loadings.  0
%   corresponds to using only initial loadings with no update
%
% 'UpsTol' (1e-6) - stop loading updates if norm(Ups(k)-Ups(k-1)) <= UpsTol
%
% 'lambda' (2.2*sqrt(n)*norminv(1-(.1/log(n))/(2*p))) - lasso penalty
%   parameter where [n,p] = size(x) in the case where observations are
%   treated as independent and p = size(x,2) and n =
%   numel(unique(clusterVar)) otherwise
%
% 'beta0' ([]) - initial coefficient estimates.  If nonempty beta0 is
%   passed, this value is used to compute initial penalty loadings and
%   'nXinit', 'demean', and 'scale' options are ignored.  beta0 should be
%   in units of passed x so that initial forecast for y is x*beta0
%
% 'iid' (0) - use iid penalty loadings (1)
%
% 'clusterVar' ([]) - cluster membership for clustered lasso
%
% IF iid = 0 AND isempty(clusterVar), HETEROSCEDASTIC LOADINGS ARE USED.
%
% 'nXinit' (ceil(log(n)) - number of variables from x to use in estimating
%   model for forming initial residuals.  The variables selected are the
%   nXinit most highly correlated to the outcome.
%
% 'demean' (0) - demean y and x before applying lasso.  The lasso and
%   post-lasso intercept estimates are returned as alasso and apost when
%   demean = 1 is used.  
%
% 'standard' (0) - standardize y and x before applying lasso
%   estimation.  Returned estimates are in terms of original units.  The
%   way this is written now, standard = 1 will break if a column of ones
%   or any other vaiable with no variation is included in x.  The lasso and
%   post-lasso intercept estimates are returned as alasso and apost when
%   standard = 1 is used.  Standardization is with respect to estimate of 
%   standard deviation under independence across observations without 
%   regard to clusterVar
%
% 'nonpen' ([]) - vector with column indexes of x that are not to be
%   penalized during estimation.  Note that this vector is not used during
%   the initial penalty loading formation in the current formulation.  This
%   option will be very inefficient for dealing with a large number of
%   unpenalized coefficients (e.g. fixed effects)
%
% 'partial' (ones(n,1)) - n x k matrix containing variables to partial out
%   ex ante.  
%   IF 'PARTIAL' IS NOT EMPTY, 'STANDARD' AND 'DEMEAN' ARE BOTH SET TO 0.  
%   ONLY THE VARIABLES IN 'PARTIAL' ARE PARTIALED OUT.  PARTIAL SHOULD THUS
%   INCLUDE AN INTERCEPT ALMOST ALWAYS.  PARTIAL SHOULD NOT OVERLAP WITH x.  
%   COEFFICIENTS FOR COLUMNS IN 'PARTIAL' ARE RETURNED IN bpart.
%   Maybe should adjust p to reflect dimension of 'partial' - especially
%   when this dimension is large.  This adjustment is not implemented.
%
% 'lassoresid' (1) - use lasso residuals in forming penalty loadings (1) or
%   use post-lasso residuals in forming penalty loadings (0)

%% Process inputs
[n,p] = size(x);  % Size of x
[Verbose,MaxIter,UpsTol,lambda,beta0,iid,clusterVar,nXinit,demean,standard,nonpen,partial,lassoresid] ...
    = process_options(varargin,'Verbose',0,'MaxIter',15,'UpsTol',1e-6,...
    'lambda',[],'beta0',[],'iid',0,'clusterVar',[],'nXinit',[],'demean',0,...
    'standard',0,'nonpen',[],'partial',ones(n,1),'lassoresid',1);

if isempty(clusterVar)
    neff = n;  % effective sample size
else
    neff = numel(unique(clusterVar));  % effective sample size
end

if isempty(nXinit)
    nXinit = min(p,ceil(log(neff)));
end
if isempty(lambda)
    lambda = 2.2*sqrt(neff)*norminv(1-(.1/log(neff))/(2*p));
end

%% Standardize data
if isempty(partial)
    if standard == 1
        [y,my,sy] = zscore(y);  % standardize y and return mean and standard deviation
        [x,mx,sx] = zscore(x);  % standardize x and return 1 x p mean and standard deviation vectors
    end
    if standard == 0 && demean == 1
        iota = ones(n,1);  % column of ones
        my = mean(y);  % scalar, mean of y
        mx = mean(x);  % 1 x p vector, mean of x
        y = y - iota*(iota\y);  % demean outcome
        x = x - iota*(iota\x);  % demean design
    end
else
    deltaY = partial\y;  % k x 1 vector of partial on y coefficients
    deltaX = partial\x;  % k x p vector of partial on x coefficients
    y = y - partial*deltaY;  % Partial out "partial" from y
    x = x - partial*deltaX;  % Partial out "partial" from x
end

%% Form initial residuals
if isempty(beta0)
    if nXinit > 0
        cyx = corr(x,y);  % p x 1 vector of correlations between x and y
        cyx = abs(cyx);  % magnitude of correlation
        sortcyx = sort(cyx,'descend');  % sort correlations
        cutcyx = sortcyx(nXinit);  % order statistic
        useinit = cyx >= cutcyx;  % columns to use for forming initial residuals
        e = y-x(:,useinit)*(x(:,useinit)\y);  % initial residual estimates
    else
        e = y;  % initial residual estimates
    end
else
    e = y - x*beta0;  % initial residual estimates
end

%% Independent cases
if isempty(clusterVar)
%% Main lasso step for iid case
    if iid == 1
        Ups0 = std(e)*std(x);  % Penalty loadings
        Ups0(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
        blasso = lassoCDFu(x, y , lambda, Ups0, 'Verbose', Verbose, 'beta', beta0);  % Initial lasso estimates
        if lassoresid
            e = y - x*blasso;   % lasso residuals
        else
            use = abs(blasso) > 0;
            e = y-x(:,use)*(x(:,use)\y); % post-lasso residuals
        end
        
        % loading updates
        kk = 1;
        Ups1 = std(e)*std(x);  % new loadings
        Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
        while kk <= MaxIter
            d0 = norm(Ups0-Ups1);  % how close are loadings from step kk-1 and those that will be used at step kk
            blasso = lassoCDFu(x, y , lambda, Ups1 , 'Verbose', Verbose, 'beta', blasso);
            if lassoresid
                e = y - x*blasso;  % lasso residuals
                Ups0 = Ups1;   % Replace old loadings
                Ups1 = std(e)*std(x); % new loadings
                Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
                kk = kk+1;  % increase counter
                d1 = norm(Ups0 - Ups1);  % how close are loadings from step kk and those that will be used at step kk+1
                if d1 == d0   % Can have a case where loadings aren't close to each other
                              % but keep cycling back and forth from adding and subtracting variables.  
                              % This clause should prevent that.
                    Ups0 = Ups1;
                end
            else
                use = abs(blasso) > 0;
                if sum(use) < n-1
                    e = y-x(:,use)*(x(:,use)\y);  % post-lasso residuals
                    Ups0 = Ups1;  % replace old loadings
                    Ups1 = std(e)*std(x);  % new loadings
                    Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
                    kk = kk+1;  % increase counter
                    d1 = norm(Ups0 - Ups1);  % how close are loadings from step kk and those that will be used at step kk+1
                    if d1 == d0   % Can have a case where loadings aren't close to each other
                                  % but keep cycling back and forth from adding and subtracting variables.  
                                  % This clause should prevent that.
                        Ups0 = Ups1;
                    end
                else
                    kk = MaxIter + 1;  % if you've got as many variables 
                      % as observations - assuming there's an intercept 
                      % that's already been partialed out -, stop because
                      % you're doing post-lasso here.
                end
            end
            if d0 < UpsTol  % If Ups_{kk-1} and Ups_{kk} were close, stop.
                kk = MaxIter + 1;
            end            
        end    
%% Main lasso step for inid case
    else        
        Syx = x.*(e*ones(1,p));  % Make score matrix.  This will be slow when n and p are big.  Could speed up with some switches.
        Ups0 = std(Syx);  % Penalty loadings
        Ups0(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
        blasso = lassoCDFu(x, y , lambda, Ups0, 'Verbose', Verbose, 'beta', beta0);  % Initial lasso estimates
        if lassoresid
            e = y - x*blasso;   % lasso residuals
        else
            use = abs(blasso) > 0;
            e = y-x(:,use)*(x(:,use)\y); % post-lasso residuals
        end
        
        % loading updates
        kk = 1;
        Syx = x.*(e*ones(1,p));
        Ups1 = std(Syx);  % new loadings
        Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
        while kk <= MaxIter
            d0 = norm(Ups0-Ups1);  % how close are loadings from step kk-1 and those that will be used at step kk
            blasso = lassoCDFu(x, y , lambda, Ups1 , 'Verbose', Verbose, 'beta', blasso);
            if lassoresid
                e = y - x*blasso;  % lasso residuals
                Ups0 = Ups1;   % Replace old loadings
                Syx = x.*(e*ones(1,p));  
                Ups1 = std(Syx); % new loadings
                Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
                kk = kk+1;  % increase counter
                d1 = norm(Ups0 - Ups1);  % how close are loadings from step kk and those that will be used at step kk+1
                if d1 == d0   % Can have a case where loadings aren't close to each other
                              % but keep cycling back and forth from adding and subtracting variables.  
                              % This clause should prevent that.
                    Ups0 = Ups1;
                end
            else
                use = abs(blasso) > 0;
                if sum(use) < n-1
                    e = y-x(:,use)*(x(:,use)\y);  % post-lasso residuals
                    Ups0 = Ups1;  % replace old loadings
                    Syx = x.*(e*ones(1,p));  
                    Ups1 = std(Syx);  % new loadings
                    Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
                    kk = kk+1;  % increase counter
                    d1 = norm(Ups0 - Ups1);  % how close are loadings from step kk and those that will be used at step kk+1
                    if d1 == d0   % Can have a case where loadings aren't close to each other
                                  % but keep cycling back and forth from adding and subtracting variables.  
                                  % This clause should prevent that.
                        Ups0 = Ups1;
                    end
                else
                    kk = MaxIter + 1;  % if you've got as many variables 
                      % as observations - assuming there's an intercept 
                      % that's already been partialed out -, stop because
                      % you're doing post-lasso here.
                end
            end
            if d0 < UpsTol  % If Ups_{kk-1} and Ups_{kk} were close, stop.
                kk = MaxIter + 1;
            end            
        end
    end
    
%% Clustered case
else
    Dt = dummyvar(categorical(clusterVar));  % Create cluster membership dummies.  This will be slow when n and number of clusters if big.  Could speed up with some switches.
    Syx = x.*(e*ones(1,p));  % Make score matrix.  This will be slow when n and p are big.  Could speed up with some switches.
    Ups0 = std(Dt'*Syx)';  % Penalty loadings
    Ups0(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
    blasso = lassoCDFu(x, y , lambda, Ups0, 'Verbose', Verbose, 'beta', beta0);  % Initial lasso estimates
    if lassoresid
        e = y - x*blasso;   % lasso residuals
    else
        use = abs(blasso) > 0;
        e = y-x(:,use)*(x(:,use)\y); % post-lasso residuals
    end
    
    % loading updates
    kk = 1;
    Syx = x.*(e*ones(1,p));
    Ups1 = std(Dt'*Syx)';  % new loadings
    Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
    while kk <= MaxIter && norm(Ups0-Ups1)>UpsTol  %norm condition missing
        d0 = norm(Ups0-Ups1);  % how close are loadings from step kk-1 and those that will be used at step kk
        blasso = lassoCDFu(x, y , lambda, Ups1 , 'Verbose', Verbose, 'beta', blasso);
        if lassoresid
            e = y - x*blasso;  % lasso residuals
            Ups0 = Ups1;   % Replace old loadings
            Syx = x.*(e*ones(1,p));
            Ups1 = std(Dt'*Syx)';  % new loadings
            Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize            
            kk = kk+1;  % increase counter
            d1 = norm(Ups0 - Ups1);  % how close are loadings from step kk and those that will be used at step kk+1
            if d1 == d0   % Can have a case where loadings aren't close to each other
                % but keep cycling back and forth from adding and subtracting variables.
                % This clause should prevent that.
                Ups0 = Ups1;
            end
        else
            use = abs(blasso) > 0;
            if sum(use) < n-1
                e = y-x(:,use)*(x(:,use)\y);  % post-lasso residuals
                Ups0 = Ups1;  % replace old loadings
                Syx = x.*(e*ones(1,p));
                Ups1 = std(Dt'*Syx)';  % new loadings
                Ups1(nonpen) = 0;  % Kill penalization on any variables that you don't want to penalize
                kk = kk+1;  % increase counter
                d1 = norm(Ups0 - Ups1);  % how close are loadings from step kk and those that will be used at step kk+1
                if d1 == d0   % Can have a case where loadings aren't close to each other
                    % but keep cycling back and forth from adding and subtracting variables.
                    % This clause should prevent that.
                    Ups0 = Ups1;
                end
            else
                kk = MaxIter + 1;  % if you've got as many variables
                % as observations - assuming there's an intercept
                % that's already been partialed out -, stop because
                % you're doing post-lasso here.
            end
        end
        if d0 < UpsTol  % If Ups_{kk-1} and Ups_{kk} were close, stop.
            kk = MaxIter + 1;
        end
    end
end
use = abs(blasso) > 0;
Ups = Ups0;

bpost = zeros(p,1);
bpost(use) = x(:,use)\y;

alasso = [];
apost = [];
bpartlasso = [];
bpartpost = [];

if isempty(partial)
    if standard == 1
        blasso = sy*(blasso./(sx'));
        bpost = sy*(bpost./(sx'));
        alasso = my - mx*blasso;
        apost = my - mx*bpost;
    end
    if standard == 0 && demean == 1
        alasso = my - mx*blasso;
        apost = my - mx*bpost;
    end
else
    bpartlasso = deltaY - deltaX*blasso;
    bpartpost = deltaY - deltaX*bpost;
end

end   

