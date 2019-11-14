function [B, B0, lamgrid] = rankflassopath(y,X,lambda2,L,eps,printitn)
% b = rankflassopath(y,X,lambda2,...)
% Computes the rank fused-Lasso regression estimates for given fused
% penalty value lambda_2 and for a range of lambda_1 values
%
% INPUT: 
%   y       : numeric response N x 1 vector (real/complex)
%   X       : numeric feature  N x p matrix (real/complex)
%   lambda2 : positive penalty parameter for the fused Lasso penalty term
%   L       : number of grid points for lambda1 (Lasso penalty)
%   eps     : Positive scalar, the ratio of the smallest to the 
%             largest Lambda value in the grid. Default is eps = 10^-4. 
% printitn : print iteration number (default = 0, no printing)
% OUTPUT:
%       B: Fitted rank fused-Lasso regression coefficients, a p-by-(L+1) matrix, 
%           where p is the number of predictors (columns) in X, and L is 
%           the  number of Lambda values.
%       B0: estimates values of intercepts
%       lamgrid: = lambda parameters 
%
% version: Sep 2, 2018 
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, p] = size(X);

intcpt = false; 

if nargin <  6
    printitn = 0;
end
if nargin < 5 || isempty(eps)
    eps  = 10^-3; 
end
if nargin < 4 || isempty(L)
    L = 120;
end

B=repmat((1:n)',1,n);
A=B';
a=A(A<B);
b=B(A<B);

D=-eye(p-1);
D(p:p:(p-1)^2)=1; 
onev = [zeros(p-2,1);1];
D = [D onev]; 
ytilde = [y(a)-y(b); zeros(p-1,1)];
Xtilde = [X(a,:)- X(b,:);lambda2*D];

lam0= norm(Xtilde'*sign(ytilde),'inf');

lamgrid = eps.^((0:1:L)/L)*lam0;
B = zeros(p,L+1); 
B0 = zeros(1,L+1);
binit = zeros(p,1);

if printitn > 0
    fprintf('rankflassopath: starting iterations\n');
end

for jj=1:(L+1) 
   B(:,jj) =ladlasso(ytilde,Xtilde,lamgrid(jj),intcpt,binit,printitn); 
   binit = B(:,jj); 
   if printitn >0 
        fprintf(' . ');
   end
end

B(abs(B)<1e-7)=0;


