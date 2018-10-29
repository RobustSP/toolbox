function [B, B0, stats] = ranklassopath(y, X, L, eps,reltol,printitn)
% [B, B0, stats] = ranklassopath(y, X,...)
% ranklassopath computes the rank LAD-Lasso regularization path (over grid 
% of penalty parameter values). Uses IRWLS algorithm.  
% INPUT: 
%        y: Numeric data vector of size N x 1 (output, respones)
%        X: Numeric data matrix of size N x p. Each row represents one 
%           observation, and each column represents one predictor (feature). 
%        L: Positive integer, the number of lambda values on the grid to be  
%           used. The default is L=120. 
%      eps: Positive scalar, the ratio of the smallest to the 
%           largest Lambda value in the grid. Default is eps = 10^-3. 
%   reltol: Convergence threshold for IRWLS. Terminate when successive 
%          estimates differ in L2 norm by a rel. amount less than reltol.
%  printitn: print iteration number (default = 0, no printing)
% OUTPUT:
%        B: Fitted RLAD-Lasso regression coefficients, a p-by-(L+1) matrix, 
%           where p is the number of predictors (columns) in X, and L is 
%           the  number of Lambda values.
%       B0: estimates values of intercepts
%    stats: structure with following fields: 
%           Lambda = lambda parameters in ascending order
%           GMeAD = Mean Absolute Deviation (MeAD) of the residuals
%           gBIC = generalized Bayesian information criterion (gBIC) value  
%                for each lambda parameter on the grid. 
%
% version: Aug 31, 2018
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, p] = size(X);
intcpt = false; 
if nargin <  6
    printitn = 0;
end

if nargin <  5 || isempty(reltol)
% sometimes IRWLS gets stuck at non-optimal points. Setting
% relative tolerance low, yields sometimes better performance  
    reltol = 1.0e-7;
end


if nargin < 4 || isempty(eps) 
    eps  = 10^-3; 
end

if  ~isscalar(eps) || ~isreal(eps) || ~isfinite(eps) || eps < 0 
    error('rladlassopath: invalid value of eps');
end

if nargin < 3 || isempty(L) 
    L = 120;
end

if  ~isscalar(L) || ~isreal(L) || ~isfinite(L) || L < 0 
    error('rladlassopath: invalid value of L');
end

B=repmat((1:n)',1,n);
A=B';
a=A(A<B);
b=B(A<B);
Xtilde = X(a,:)- X(b,:); 
ytilde = y(a)-y(b);
lam0= norm(Xtilde'*sign(ytilde),'inf');

lamgrid = eps.^((0:1:L)/L)*lam0;
B = zeros(p,L+1); 
B0 = zeros(1,L+1);
binit = zeros(p,1);

for jj=1:(L+1) 
   B(:,jj)= ladlasso(ytilde,Xtilde,lamgrid(jj),intcpt,binit,reltol,printitn); 
   binit = B(:,jj); 
   r = y-X*binit;   
   if isreal(X) 
        B0(jj) = median((r(a)+r(b))/2);
   else
        B0(jj) = spatmed((r(a)+r(b))/2);
   end
end
B(abs(B)<1e-7)=0;
stats.DF = sum(abs(B)~=0); 

%% Compute the generalized BIC (gBIC) values
Rmat   = repmat(ytilde,1,L+1)-Xtilde*B;     % matrix of residuals 
N=n*(n-1)/2;
stats.GMeAD = (sqrt(pi)/2)*mean(abs(Rmat)); % Gini's dispersion  
stats.GMeAD = stats.GMeAD.*sqrt(n./(n-stats.DF-1));
stats.gBIC = 2*n*log(stats.GMeAD) + stats.DF.*log(n);
stats.Lambda = lamgrid;

