function  [B, B0, stats] = hublassopath(y,X,c,L,eps,intcpt,reltol,printitn)
% [B, B0, stats] = hublassopath(y, X,...)
% hublassopath computes the M-Lasso regularization path (over grid 
% of penalty parameter values) using Huber's loss function  
% INPUT: 
%        y: Numeric data vector of size N x 1 (output, respones)
%        X: Numeric data matrix of size N x p. Each row represents one 
%           observation, and each column represents one predictor (feature)
%           columns are standardized to unit length.
%        c: Threshold constant of Huber's loss function (optional;
%            otherwise use defaul falue)
%   intcpt: Logical (true/false) flag to indicate if intercept is in the 
%           regression mode. Default is true.
%      eps: Positive scalar, the ratio of the smallest to the 
%           largest Lambda value in the grid. Default is eps = 10^-3. 
%       L : Positive integer, the number of lambda values EN/Lasso uses.  
%           Default is L=120. 
%   reltol : Convergence threshold for IRWLS. Terminate when successive 
%           estimates differ in L2 norm by a rel. amount less than reltol.
% printitn: print iteration number (default = 0, no printing)
% OUTPUT:
%   B    : Fitted M-Lasso regression coefficients, a p-by-(L+1) matrix, 
%          where p is the number of predictors (columns) in X, and L is 
%          the  number of Lambda values.
%     B0 : estimates values of intercepts
%  stats  : structure with following fields: 
%           Lambda = lambda parameters in ascending order
%           sigma = estimates of the scale (a (L+1) x 1 vector)
%           gBIC = generalized Bayesian information criterion (gBIC) value  
%                 for each lambda parameter on the grid. 
% version: Sep 20, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,p] = size(X);

realdata = true;
if ~isreal(y) 
    realdata=false; 
end

if nargin <  8
    printitn = 0;
end

if nargin <  7 || isempty(reltol)
% sometimes IRWLS gets stuck at non-optimal points and therefore setting
% relative tolerance low, yields sometimes better performance  
    reltol = 1.0e-5;
end

if nargin < 6  || isempty(intcpt)
    intcpt=true; 
end 

if nargin < 5 || isempty(eps)
    eps  = 10^-3; 
end

if nargin < 4 || isempty(L)
    L = 120;
end

if nargin < 3 || isempty(c)
   if realdata
       c = 1.3415;
   else
       c = 1.215;
   end % Default: approx 95 efficiency for Gaussian errors
end

[locy, sig0] = hubreg(y,ones(n,1),c);

if intcpt % if intercept is in the model, center the data
    ync = y; Xnc = X;
    meanX = mean(X);
    X = bsxfun(@minus, X, meanX);
    y = y-locy;
end         

% standardize the predictors to unit norm columns
sdX = sqrt(sum(X.*conj(X)));
X = bsxfun(@rdivide, X, sdX);

% compute the smallest penalty value yielding a zero solution
yc = psihub(y/sig0,c)*sig0;
lam0 =norm(X'*yc,'inf');

lamgrid = eps.^((0:1:L)/L)*lam0;
B = zeros(p,L+1); 
sig = zeros(1,L+1); sig(1) = sig0;
for jj=1:L
   [B(:,jj+1), sig(jj+1)] = hublasso(y,X,c,lamgrid(jj+1),B(:,jj),sig(jj),reltol,printitn); 
end
B(abs(B)<5e-8)=0;
DF = sum(abs(B)~=0);
con = sqrt((n./(n-DF-1)));

if n > p 
    stats.gBIC = 2*n*log(sig.*con) + DF.*log(n); 
else 
    stats.gBIC = [];
end

B = bsxfun(@rdivide,B,sdX'); 
% compute the intercept if in the model 
if intcpt 
    B0  = locy - meanX*B;
else
    B0 = [];
end

stats.sigma = sig;
stats.Lambda = lamgrid;

