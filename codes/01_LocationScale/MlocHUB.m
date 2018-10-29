function [mu_hat] = MlocHUB(y,c)
% Mloc_HUB computes Huber's M-estimate of
% location, i.e.,
%
% mu_hat = arg min_mu SUM_i rho_HUB(y_i - mu)
%
%
%   INPUTS: 
%           y: real valued data vector of size N x 1
%           c: tuning constant c>=0 
%
%   OUTPUT:  
%           mu_hat: Huber's M-estimate of location
% version: Sep 25, 2018
% authors: Michael Muma 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
   c = 1.345; % default tuning for 95 percent efficiency under the Gaussian model      
end

% previously computed scale estimate
    sigma_0 = madn(y); 

% initial robust location estimate 
    mu_n = median(y);

% maximum number of iterations
    max_iters = 1000;

% convergence error tolerance
    tol_err = 1.0e-5;
    
% iteration counter    
    n = 0;

while  n<=max_iters
    
 w_n =  whub(abs(y-mu_n)/sigma_0,c); % compute weights
 mu_n_plus1 = sum(w_n.*y)/(sum(w_n)); % compute weighted average
 
 if abs(mu_n_plus1-mu_n)/sigma_0 > tol_err % breaking condition
        mu_n = mu_n_plus1; % update estimate of mean
        n = n+1; % increment iteration counter      
 else
     break
 end
end

mu_hat = mu_n; % final estimate

end
