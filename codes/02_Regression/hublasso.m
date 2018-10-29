function  [b0, sig0, psires] = hublasso(y,X,c,lambda,b0,sig0,reltol,printitn)
% [B, stats] = hublasso(y, X,c,lambda,b0,sig0,...)
% hublasso computes the M-Lasso estimate for a given penalty parameter 
% using Huber's loss function 
% INPUT: 
%        y: Numeric data vector of size N x 1 (output,respones)
%        X: Numeric data matrix of size N x p (inputs,predictors,features). 
%           Each row represents one observation, and each column represents 
%           one predictor
%        c: Threshold constant of Huber's loss function 
%   lambda: positive penalty parameter value
%       b0: numeric initial start of the regression vector
%     sig0: numeric positive scalar, initial scale estimate.  
%   reltol: Convergence threshold. Terminate when successive 
%           estimates differ in L2 norm by a rel. amount less than reltol.
%           Default is 1.0e-5
% printitn: print iteration number (default = 0, no printing)
% OUTPUT:
%       b0: regression coefficient vector estimate
%     sig0: estimate of the scale 
%   psires: pseudoresiduals
%
% version: Aug 31, 2018
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,p] = size(X);

realdata = true;
if ~isreal(y) 
    realdata=false; 
end

if nargin <  7
    printitn = 0;
end

if nargin <  6 || isempty(reltol)
% sometimes IRWLS gets stuck at non-optimal points and therefore setting
% relative tolerance low, yields sometimes better performance  
    reltol = 1.0e-5;
end

if nargin < 3 || isempty(c) 
    if  realdata
        c = 1.345; 
    else
        c = 1.1214; 
    end
end

csq = c^2;
if realdata % real case 
    qn = chi2cdf(csq,1);
    alpha = chi2cdf(csq,3)+csq*(1-qn); % consistency factor for scale 
else % complex case
    qn = chi2cdf(2*csq,2);
    alpha = chi2cdf(2*csq,4)+csq*(1-qn); % consistency factor for scale 
end
iterMAX = 500;
con  = sqrt(n*alpha);
betaold = b0;
normb0 = norm(b0);

for iter = 1:iterMAX
 
    r = y - X*b0;
    psires = psihub(r/sig0,c)*sig0;
    sig1 = norm(psires)/con;  % update the scale
    
    crit2 = abs(sig0-sig1); 
    for jj=1:p           
        psires = psihub(r/sig1,c)*sig1;  % Update the pseudo-residual             
        b0(jj) = SoftThresh(b0(jj)+X(:,jj)'*psires,lambda);
        r = r+X(:,jj)*(betaold(jj)-b0(jj));
    end
    
    normb = norm(b0);
    crit = sqrt(normb0^2 + normb^2 - 2*real(betaold'*b0))/normb;

    if printitn > 0 
        r = (y - X*b0)/sig1;
        objnew = (sig1)*sum(rhofun(r,c)) + (beta/2)*n*sig1 + lambda*sum(abs(b0));
        fprintf('Iter %3d obj =  %10.5f crit = %10.5f crit2 = %10.5f\n', ... 
            iter,objnew,crit,crit2);
    end

     if crit < reltol
        break
    end
    
    sig0 = sig1;
    betaold = b0;
    normb0 = normb;   
end

if printitn > 0
% Check the M-Lasso estimating equations
    b0(abs(b0)<5e-9)=0;
    r = y - X*b0;
    s = sign(b0);
    ind  = 1:p;
    ind2 = ind(s==0);
    psires = psihub(r/sig1,c)*sig1;  
    s(ind2)=  X(:,ind2)'*psires/lambda;
    fpeq = -X'*psires + lambda*s;       %  FP equation equal to zero
    fprintf('lam = %.4f it = %d norm(FPeq1)= %.12f abs(FPeq2)=%.12f\n', ...
        lambda,iter, norm(fpeq),sig1-norm(psires)/con);
end

