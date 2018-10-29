function [b1,iter] = rladreg(y,X,b0,printitn)
% [b1, iter] = rladreg(y,X,...)
% ladreg computes the LAD regression estimate 
% INPUT:
%         y: numeric response N x 1 vector (real/complex)
%         X: numeric feature  N x p matrix (real/complex)
%        b0: numeric optional initial start of the regression vector for 
%            IRWLS algorithm. If not given, we use LSE (when p>1).
%  printitn: print iteration number (default = 0, no printing) and
%            other details 
% OUTPUT:
%        b1: numeric the regression coefficient vector
%      iter: (numeric) # of iterations (given when IRWLS algorithm is used)
%
% version  : Sep 7, 2018 (Dependencies: uses ladlasso function
% authors  : Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, p] = size(X);


if nargin < 4
    printitn = 0;
end


if nargin<3 
    b0  = [ones(n,1) X] \ y;
    b0  = b0(2:end);
end

B=repmat((1:n)',1,n);
A=B';
a=A(A<B);
b=B(A<B);
Xtilde = X(a,:)- X(b,:); 
ytilde = y(a)-y(b);

if p==1
     b1 = wmed(ytilde./Xtilde,abs(Xtilde));
     iter = [];
else
     [b1, iter] = ladlasso(ytilde,Xtilde,0,b0,printitn);
end




