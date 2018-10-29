function [C,smed0] = signcm(X,center);
% C = scm(X,center);
% calculates the spatial sign covariance matrix (SCM). 
%
% INPUT 
%         X: Numeric data matrix of size N x p. Each row represents one 
%           observation, and each column represents one variable.
%    center: logical (true/false). If true, then center the data using
%            spatial median. Default is false
% OUTPUT: 
%         C: spatial sign covariance matrix
%     smed0: spatial median (computed only if center = true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


EPS = 1.0e-6;

if nargin < 2 
    center = false;
    smed0 = [];
end

if center 
    smed0 = spatmed(X);
    X = bsxfun(@minus,X,smed0);
end

len = sqrt(sum(X.*conj(X),2));

X(len~=0,:)= X;
len = len(len~=0);
n = size(X,1);
len(len<EPS)= EPS;
X = bsxfun(@rdivide, X, len); 
C = X.'*conj(X)/n;


