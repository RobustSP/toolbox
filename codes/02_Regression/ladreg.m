function [b1, iter] = ladreg(y,X,intcpt,b0,printitn)
% [b1, iter] = ladreg(y,X,intcpt,...)
% ladreg computes the LAD regression estimate 
% INPUT: 
%         y: numeric response N x 1 vector (real/complex)
%         X: numeric feature  N x p matrix (real/complex)
%    intcpt: (logical) flag to indicate if intercept is in the model
%        b0: numeric optional initial start of the regression vector for 
%            IRWLS algorithm. If not given, we use LSE (when p>1).
%  printitn: print iteration number (default = 0, no printing) and
%            other details 
% OUTPUT:
%        b1: (numberic) the regression coefficient vector
%      iter: (numeric) # of iterations (given when IRWLS algorithm is used)
% 
% version: Aug 31, 2018 
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    printitn = 0;
end

if nargin < 4 
    b0 = [];
end

if nargin < 3
    error('ladgreg: not enough input arguments');
end

[b1,iter] = ladlasso(y,X,0,intcpt,b0,printitn);



