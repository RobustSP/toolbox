function [b1, sig1, iter] = hubreg(y,X,c,sig0,b0,printitn)
% [b1, sig1, iter] = hubreg(y,X,...)
% hubreg computes the joint M-estimates of regression and scale using 
% Huber's criterion. Function worsk for both real- and complex-valued data.
%
% INPUT: 
%        y: Numeric data vector of size N x 1 (output, respones)
%        X: Numeric data matrix of size N x p. Each row represents one 
%           observation, and each column represents one predictor (feature). 
%           If the model has an intercept, then first column needs to be a  
%           vector of ones. 
%         c: numeric threshold constant of Huber's function
%      sig0: (numeric) initial estimator of scale [default: SQRT(1/(n-p)*RSS)]
%        b0: initial estimator of regression (default: LSE)  
%  printitn: print iteration number (default = 0, no printing)
% OUTPUT:
%       b1: the regression coefficient vector estimate 
%     sig1: the estimate of scale 
%     iter: the # of iterations 
%
% version: Sep 2, 2018
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, p] = size(X);

realdata = true;
if ~isreal(y) 
    realdata=false; 
end

if nargin < 6 
    printitn = 0;
end

if nargin < 3 || isempty(c)
   if realdata
       c = 1.3415;
   else
       c = 1.215;
   end % Default: approx 95 efficiency for Gaussian errors
end

if nargin < 5 || isempty(b0)
% if initial start not given, use the LSE    
    b0 = X \ y;
end

if nargin < 4 || isempty(sig0)
% initial estimate of scale -->  (unbiased) RSS    
    sig0 =  norm(y-X*b0)/sqrt(n-p);
end

csq = c^2; 
if realdata % real case 
    qn = chi2cdf(csq,1);
    alpha = chi2cdf(csq,3)+csq*(1-qn); % consistency factor for scale 
else % complex case
    qn = chi2cdf(2*csq,2);
    alpha = chi2cdf(2*csq,4)+csq*(1-qn); % consistency factor for scale 
end
  
ITERMAX = 2000;
ERRORTOL = 1e-5; % ERROR TOLERANCE FOR HALTING CRITERION

Z = pinv(X);
con  = sqrt((n-p)*alpha);

for iter=1:ITERMAX
   
   % STEP 1: update residual  
   r = y - X*b0;
   psires = psihub(r/sig0,c)*sig0;
   
   % STEP 2 Update the scale   
   sig1= norm(psires)/con;
      
   % STEP 3 Update the pseudo-residual 
   psires = psihub(r/sig1,c)*sig1;  
   
   % STEP 4 regresses X on pseudoresidual
   update = Z*psires;
  
   % STEP 5 update the beta
   b1 = b0 + update; 
  
   % STEP 6 (Check convergence)    
   crit2 = norm(update)/norm(b0);
  
   if mod(iter,printitn)==0
        fprintf('hubreg: crit(%4d) = %.9f\n',iter,crit2); 
   end
   
   if crit2 < ERRORTOL
      break 
   end   
   
   b0 = b1; sig0 = sig1; 
end

if iter == ITERMAX
    fprintf('error!!! MAXiter = %d crit2 = %.7f\n',iter,crit2)
end


