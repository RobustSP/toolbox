function [B, stats] = ladlassopath(y,X,L,eps,intcpt,reltol,printitn)
% [B, stats] = ladlassopath(y, X,...)
% ladlassopath computes the LAD-Lasso regularization path (over grid 
% of penalty parameter values). Uses IRWLS algorithm.  
% INPUT: 
%       y : Numeric data vector of size N x 1 (output, respones)
%       X : Numeric data matrix of size N x p. Each row represents one 
%           observation, and each column represents one predictor (feature). 
%   intcpt: Logical (true/false) flag to indicate if intercept is in the 
%           regression model
%      eps: Positive scalar, the ratio of the smallest to the 
%           largest Lambda value in the grid. Default is eps = 10^-3. 
%       L : Positive integer, the number of lambda values EN/Lasso uses.  
%           Default is L=120. 
%   reltol : Convergence threshold for IRWLS. Terminate when successive 
%           estimates differ in L2 norm by a rel. amount less than reltol.
% printitn: print iteration number (default = 0, no printing)
% OUTPUT:
%   B    : Fitted LAD-Lasso regression coefficients, a p-by-(L+1) matrix, 
%          where p is the number of predictors (columns) in X, and L is 
%          the  number of Lambda values. If intercept is in the model, then
%          B is (p+1)-by-(L+1) matrix, with first element the intercept.
%  stats  : structure with following fields: 
%           Lambda = lambda parameters in ascending order
%           MeAD = Mean Absolute Deviation (MeAD) of the residuals
%           gBIC = generalized Bayesian information criterion (gBIC) value  
%                 for each lambda parameter on the grid. 
% version: Aug 31, 2018
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[n, p] = size(X);

if nargin <  7
    printitn = 0;
end

if nargin <  6 || isempty(reltol)
% sometimes IRWLS gets stuck at non-optimal points and therefore setting
% relative tolerance low, yields sometimes better performance  
    reltol = 1.0e-6;
end

if nargin < 5  || isempty(intcpt)
    intcpt=true; 
end 

if ~islogical(intcpt)
   error('ladlassopath: input ''intcpt'' needs to be ''false'' or ''true''');
end

if nargin < 4 || isempty(eps)
    eps  = 10^-3; 
end

if  ~isscalar(eps) || ~isreal(eps) || ~isfinite(eps) || eps < 0 
    error('ladlassopath: invalid value of eps');
end

if nargin < 3 || isempty(L)
    L = 120;
end

if  ~isscalar(L) || ~isreal(L) || ~isfinite(L) || L < 0 
    error('ladlassopath: invalid value of L');
end


if intcpt % if intercept is in the model
    p = p + 1;
    if isreal(y) 
        medy = median(y);
    else
        medy = spatmed(y);
    end
    yc = y - medy;
    lam0 = norm(X'*sign(yc),'inf');
else
    lam0 =norm(X'*sign(y),'inf');
end     

lamgrid = eps.^((0:1:L)/L)*lam0;
B = zeros(p,L+1); 

if intcpt 
    binit = [medy; zeros(p-1,1)]; 
else
    binit = zeros(p,1);
end
    
for jj=1:(L+1) 
   B(:,jj) = ladlasso(y,X,lamgrid(jj),intcpt,binit,reltol,printitn);
   binit = B(:,jj);
end

if intcpt
    B([false(1,L+1); abs(B(2:end,:))<1e-7])=0;
    stats.DF = sum(abs(B(2:end,:))~=0) ;
    stats.MeAD = sqrt(pi/2)*mean(abs(repmat(y,1,L+1)- [ones(n,1) X]*B)); 
    const = sqrt(n./(n-stats.DF-1));
else
    B(abs(B)<1e-7)=0;
    stats.DF = sum(abs(B)~=0);
    stats.MeAD = sqrt(pi/2)*mean(abs(repmat(y,1,L+1)- X*B)); 
    const = sqrt(n./(n-stats.DF));
end

stats.MeAD = stats.MeAD.*const;
stats.gBIC = 2*n*log(stats.MeAD) + stats.DF.*log(n); % BIC values
stats.Lambda = lamgrid;



