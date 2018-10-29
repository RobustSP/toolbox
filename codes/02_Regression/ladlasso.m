function [b1, iter] = ladlasso(y,X,lambda,intcpt,b0,reltol,printitn) 
% [b1, iter] = ladlasso(y,X,lambda,...) 
% ladlasso computes the LAD-Lasso regression estimates for given complex-  
% or real-valued data.  If number of predictors, p, is larger than one, 
% then IRWLS algorithm is used, otherwise a weighted median algorithm 
% (N > 200) or elemental fits (N<200).
% INPUT: 
%   y      : numeric response N x 1 vector (real/complex)
%   X      : numeric feature  N x p matrix (real/complex)
%   lambda : non-negative penalty parameter
%   b0     : numeric optional initial start of the regression vector for 
%           IRWLS algorithm. If not given, we use LSE (when p>1).
%   intcpt : (logical) flag to indicate if intercept is in the model
%   reltol : Convergence threshold for IRWLS. Terminate when successive 
%           estimates differ in L2 norm by a rel. amount less than reltol.
% printitn : print iteration number (default = 0, no printing)
% OUTPUT:
%   b1     : (numberic) the regression coefficient vector
%   iter   : (numeric) # of iterations (only when IRWLS algorithm is used)
%
% version: Sep 2, 2018 (Dependencies: uses ladlasso function
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, p]=size(X);

if nargin <  7
    printitn = 0;
end

if nargin <  6 || isempty(reltol)
% sometimes IRWLS gets stuck at non-optimal points. Often setting
% relative tolerance very low, yields sometimes better performance  
    reltol = 1.0e-8;
end

if nargin<4  
% Default: intercept is in the model
    intcpt=true;  
end

if ~islogical(intcpt)
   error('ladlasso: input ''intcpt'' needs to be ''false'' or ''true''');
end

if intcpt 
% add  a vector of ones as a column to X when intercept = true
   X = [ones(N,1) X];
end 

if nargin<5 || isempty(b0)
% LSE is the initial start of iteration if the initial start was not given
    b0 = X \ y;
end 

if lambda <0 || ~isreal(lambda)
      error('ladlasso: input ''lambda'' needs to be a positive scalar');
end


ITERMAX = 2000;
% we use very small error tolerance value TOL between iterations. 
iter = [];

if printitn > 0 
    fprintf('Computing the solution for lambda = %.3f\n',lambda);
end

%% the case of only one predictor 
if p==1 && ~intcpt % simple linear regression and no intercept 
        b1 = wmed([y./X; 0], [abs(X); lambda]);
elseif p==1 && isreal(y) && N < 200  && intcpt         
    if lambda==0  
        b = elemfits(X(:,2),y);            
    else
        b = elemfits([X(:,2);0],[y;lambda]);
    end
    res = sum(abs(repmat(y,1,size(b,2)) - X * b));
    [~, indx] = min(res);
    b1 = b(:,indx);
else
%%  use IRWLS always when p > 1 
    if printitn > 0
        fprintf('Starting the IRWLS algorithm..\n');
    end
    if lambda >0
        y = [y; zeros(p,1)];
        if intcpt
             X = [X; zeros(p,1) lambda*eye(p)];
        else
             X=  [X; lambda*eye(p)];
        end
    end
    
    for iter=1:ITERMAX
      
        resid = abs(y-X*b0);  
        resid(resid<.000001)=.000001;
        Xstar = bsxfun(@rdivide, X, resid);
        b1 = (Xstar'*X) \ (Xstar'*y);
 
        crit = norm(b1-b0)/norm(b0);  
        %crit = norm(b1-b0,Inf);
        if mod(iter,printitn)==0
            fprintf('ladlasso: crit(%4d) = %.9f\n',iter,crit); 
        end
        
        if crit < reltol && iter > 10
            break 
        end 
  
        b0 = b1;
    end
end

