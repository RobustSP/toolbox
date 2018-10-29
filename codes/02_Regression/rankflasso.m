function [b,iter] = rankflasso(y,X,lambda1,lambda2,b0,printitn)
% b = rankflasso(y,X,lambda1,lambda2,...)
% Computes the rank fused-Lasso regression estimates for given fused
% penalty value lambda_2 and for a range of lambda_1 values
%
% INPUT: 
%   y       : numeric response N x 1 vector (real/complex)
%   X       : numeric feature  N x p matrix (real/complex)
%   lambda1 : positive penalty parameter for the Lasso penalty term
%   lambda2 : positive penalty parameter for the fused Lasso penalty term
%   b0      : numeric optional initial start (regression vector) of 
%           iterations. If not given, we use LSE (when p>1).
% printitn : print iteration number (default = 0, no printing)
% OUTPUT:
%   b      : numeric regression coefficient vector
%   iter   : positive integer, the number of iterations of IRWLS algorithm
%
% version: Sep 2, 2018 
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, p] = size(X);

intcpt = false; 

if nargin <  6
    printitn = 0;
end

if nargin<5 || isempty(b0)
    b0  = [ones(n,1) X] \ y;
    b0  = b0(2:end);
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

%lam0= norm(Xtilde'*sign(ytilde),'inf')
if printitn > 0 
    fprintf('rankflasso: starting iterations\n');
end
[b,iter]= ladlasso(ytilde,Xtilde,lambda1,intcpt,b0,printitn); 
b(abs(b)<1e-7)=0;


