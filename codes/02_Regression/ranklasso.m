function [b1,iter] = ranklasso(y,X,lambda,b0,printitn)
% [b1, iter] = ranklasso(y,X,lambda,...)
% ladreg computes the rank (LAD-)regression estimate 
% INPUT:
%       y  : numeric data vector of size N x 1 (output, respones)
%       X  : numeric data matrix of size N x p (input, features)
%   lambda : penalty parameter (>= 0) 
%      b0  : numeric optional initial start (regression vector) of 
%           iterations. If not given, we use LSE. 
% printitn : print iteration number (default = 0, no printing) and
%            other details 
% OUTPUT:
%   b1     : numeric the regression coefficient vector
%   iter   : (numeric) # of iterations (given when IRWLS algorithm is used)
%
% version  : Sep 7, 2018 (Dependencies: uses ladlasso function
% authors  : Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(X);

intcpt = false;

if nargin < 5
    printitn = 0;
end

if nargin<4 || isempty(b0)
    b0  = [ones(n,1) X] \ y;
    b0  = b0(2:end);
end

B=repmat((1:n)',1,n);
A=B';
a=A(A<B);
b=B(A<B);
Xtilde = X(a,:)- X(b,:); 
ytilde = y(a)-y(b);

[b1,iter] = ladlasso(ytilde,Xtilde,lambda,intcpt,b0,printitn); 
