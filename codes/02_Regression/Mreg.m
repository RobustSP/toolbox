function [b1, sig] = Mreg(y,X,lossfun,b0,verbose)
% [b1, sig] = Mreg(y,X,...)
% Mreg computes the M-estimates of regression using an auxiliary scale
% estimate. It uses the iterative reweighted least squares (IRWLS) algorithm   
%
% INPUTS: 
%       y : (numeric) data vector of size N x 1 (output, response vector)
%       X : (numeric) data matrix of size N x p (input, feature matrix)
%           If the model has intercept, then first column of X should be a 
%           vector of ones. 
% lossfun : (string) either 'huber' or 'tukey' to identify the desired 
%           loss function
%      b0 : (numeric) Optional robust initial start (regression vector) of 
%           iterations. If not given, we use the LAD regression estimate 
%  verbose: (logical) true of false (default). Set as true if you wish  
%           to see convergence as iterations evolve.
% OUTPUTS:
%      
% version: Aug 31, 2018 
% Dependencies: ladreg,ladlasso, whub, and wtuk
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Compute the LAD estimate: used for estimating scale and as initial start
if nargin < 5
    verbose = false; 
end

if nargin < 4 || isempty(b0)
    b0 = ladreg(y,X,false); 
end

% Compute the auxiliary scale estimate as 
if isreal(y)  
    const =1.4815; 
else
    const =1.20112; 
end
resid = abs(y-X*b0); 
sig = const*median(resid(resid~=0)); % auxiliary scale estimate 

% Choose loss function:
if nargin == 2 || isempty(lossfun)
    lossfun = 'huber';
end

switch lossfun
    case 'huber'
        if isreal(y)
            c = 1.345; % 95 percent ARE
        else
            c =1.214; % 95 percent ARE
        end
        wfun = @(x) whub(x,c); 
    case 'tukey'
        if isreal(y)
            c = 3.4437; % 85 percent ARE
        else
            c = 3.0; % 85 percent ARE
        end
        wfun = @(x) wtuk(x,c); 
    otherwise
        error(['Mreg: input ''lossfun'' needs to be a character',...
            'string equal to ''tukey'' or ''huber''\n']);        
end

ITERMAX = 1000;
TOL = 1.0e-5;
  

if verbose
    fprintf('Mreg: iterations starting, using %s loss function \n',lossfun);
end


for iter=1:ITERMAX 
      
     resid(resid<.000001)=.000001;
     w = wfun(resid/sig);
     Xstar = bsxfun(@times, X, w);
     b1 = (Xstar'*X) \ (Xstar'*y);
     
     crit = norm(b1-b0)/norm(b0);  
     if verbose && mod(iter,1)==0
        fprintf('Mreg: crit(%4d) = %.9f\n',iter,crit); 
     end
     
     if crit < TOL
        break;
     end
    
     b0 = b1;
     resid = abs(y-X*b0);
end



