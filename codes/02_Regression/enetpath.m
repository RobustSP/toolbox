function [B, stats] = enetpath(y, X, alpha, L, eps, intcpt,printitn)
% [B, lamgrid, BIC, MSE] = enetpath(y, X,...)
% enethpath computes the elastic net (EN) regularization path (over grid 
% of penalty parameter values). Uses pathwise CCD algorithm. 
% INPUT: 
%       y : Numeric data vector of size N x 1 (output, respones)
%       X : Nnumeric data matrix of size N x p. Each row represents one 
%           observation, and each column represents one predictor (feature). 
%   intcpt: Logical flag to indicate if intercept is in the model
%  alpha  : Numeric scalar, elastic net tuning parameter in the range [0,1].
%           If not given then use alpha = 1 (Lasso)
%      eps: Positive scalar, the ratio of the smallest to the 
%           largest Lambda value in the grid. Default is eps = 10^-4. 
%       L : Positive integer, the number of lambda values EN/Lasso uses.  
%           Default is L=100. 
% printitn: print iteration number (default = 0, no printing)
% OUTPUT:
%   B    : Fitted EN/Lasso regression coefficients, a p-by-(L+1) matrix, 
%          where p is the number of predictors (columns) in X, and L is 
%          the  number of Lambda values. If intercept is in the model, then
%          B is (p+1)-by-(L+1) matrix, with first element the intercept.
%  stats  : structure with following fields: 
%           Lambda = lambda parameters in ascending order
%           MSE = Mean squared error (MSE)
%           BIC = Bayesian information criterion values for each lambda
% version: Aug 31, 2018
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, p] = size(X);

if nargin <  7
    printitn = 0;
end
if nargin < 6 || isempty(intcpt)
    intcpt = true;
end
if nargin < 5 || isempty(eps)
    eps  = 10^-3; 
end
if nargin < 4 || isempty(L)
    L = 120;
end
if nargin < 3 || isempty(alpha)
    alpha = 1;
end

if ~isscalar(alpha) || ~isreal(alpha) || ~isfinite(alpha) || ...
        alpha <= 0 || alpha > 1
    error('enetpath: invalid value of alpha\n');
end

if  ~isscalar(eps) || ~isreal(eps) || ~isfinite(eps) || eps < 0 
    error('enetpath: invalid value of eps\n');
end

if ~islogical(intcpt)
   error('enetpath: input ''intcpt'' needs to be ''false'' or ''true''');
end

if n < 3
    error('enetpath: too few observations');
end

if intcpt % if intercept is in the model, center the data
    meanX= mean(X);
    meany = mean(y);
    X = bsxfun(@minus, X, meanX);
    y = y-meany;
end     

if printitn > 0
    fprintf('enetpath: using alpha = %.1f \n', alpha);
end
sdX = sqrt(sum(X.*conj(X)));
X = bsxfun(@rdivide, X, sdX);
 
lam0 =norm(X'*y,'inf')/alpha; % smallest penalty value giving zero solution

lamgrid = eps.^((0:1:L)/L)*lam0; % grid of penalty values
B = zeros(p,L+1); 
for jj=1:L
   B(:,jj+1) = enet(y,X,B(:,jj),lamgrid(jj+1),alpha,printitn); 
end
B(abs(B)<5e-8)=0;

DF = sum(abs(B)~=0);  % non-zero values in each
if n> p 
    stats.MSE = sum(abs(repmat(y,1,L+1)- X*B).^2).*(1./(n-DF-1)); 
    stats.BIC = n*log(stats.MSE) + DF.*log(n); 
else 
    stats.MSE = []; 
    stats.BIC = [];
end

B = bsxfun(@rdivide,B,sdX'); 
if intcpt
    B= [meany - meanX*B;B];
end
stats.Lambda = lamgrid;


