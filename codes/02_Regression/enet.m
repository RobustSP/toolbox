function [beta,iter] = enet(y,X,beta,lambda,alpha,printitn)
% [beta,iter] = enet(y,X,beta,lambda,...)
% enet computes the elastic net estimator using the cyclic co-ordinate 
% descent (CCD) algorithm.
% 
% INPUT: 
%       y : (numeric) data vector of size N x 1 (output, respones)
%         if the intercept is in the model, then y needs to be centered. 
%       X : (numeric) data matrix of size N x p (input, features)
%           Columns are assumed to be standardized (i.e., norm(X(:,j))=1)
%           as well as centered (if intercept is in the model). 
%    beta : (numeric) regression vector for initial start for CCD algorithm
%  lambda : (numeric) a postive penalty parameter value 
%  alpha  : (numeric) elastic net tuning parameter in the range [0,1]. If
%           not given then use alpha = 1 (Lasso)
% printitn: print iteration number (default = 0, no printing)
% OUTPUT:
%   beta    : (numberic) the regression coefficient vector
%   iter  : (numeric) # of iterations 
% version: Aug 31, 2018
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, p] = size(X);

if nargin < 6
    printitn = 0;
end

betaold = beta;
normb0 = norm(beta);
r = y - X*beta; 

if nargin < 5 || isempty(alpha)
    alpha = 1; % Default is Lasso
end  

if ~isscalar(alpha) || ~isreal(alpha) || ~isfinite(alpha) || ...
        alpha <= 0 || alpha > 1
    error('enet: input ''alpha'' needs to be in the range [0,1]');
end

if lambda <0 || ~isreal(lambda)
      error('enet: input ''lambda'' needs to be a positive scalar');
end
  
lam1 = alpha*lambda;
lam2 = (1-alpha)*lambda;
const = 1/(1+lam2);

iterMAX = 1000;
if printitn > 0
    fprintf('enet : using penalty lambda = %.5f\n',lambda);
end

for iter = 1:iterMAX

    for jj=1:p
        beta(jj)= const*SoftThresh(beta(jj)+X(:,jj)'*r,lam1);
        r = r+X(:,jj)*(betaold(jj)-beta(jj));
    end
    
    normb = norm(beta);
    crit = sqrt( normb0^2 + normb^2 - 2*real(betaold'*beta))/normb;
    
    if mod(iter,printitn)==0 
        %penval = lambda*(alpha*sum(abs(beta))+ 0.5*(1-alpha)*norm(beta)^2);
        %objnew =  (1/2)*sum(abs(r).^2) + penval;
        fprintf('enet: %4d  crit = %.8f\n',iter,crit);
        %fprintf('enet: %4d  obj =  %.14f crit = %.8f\n',iter,objnew, crit);
    end

    if crit < 1.0e-4
        break
    end
    
    betaold = beta;
    normb0 = normb;
end

