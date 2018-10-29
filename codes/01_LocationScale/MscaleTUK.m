function sigma_hat = MscaleTUK(y,c)
% Mscale_TUK computes Tukey's M-estimate of
% scale.
%
%
%   INPUTS: 
%           y: real valued data vector of size N x 1
%           c: tuning constant c>=0 
%
%   OUTPUT:  
%           sigma_hat: Tukey's M-estimate of scale
% version: Sep 26, 2018
% authors: Michael Muma 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2
    c = 4.685; % default tuning for 95 percent efficiency under the Gaussian model      
end


% initial scale estimate
    sigma_n = madn(y);

% subtract previously computed location
    mu_hat = median(y);
    y = y-mu_hat;     

% maximum number of iterations
    max_iters = 1000;

% convergence error tolerance
    tol_err = 1e-5;
    
% iteration counter    
    n = 0;

% length of sample    
    N = length(y);

% consistency with the standard deviation at the Gaussian
    rng(1);
    u = randn(10000,1);
    delta = mean(rhotuk(u,c));


while  n<=max_iters   
    
w_n = wtuk(abs(y)/sigma_n,c);

sigma_n_plus1 = sqrt(1/(N*delta)*sum(w_n.*y.^2)); 

if abs(sigma_n_plus1/sigma_n-1)>tol_err
sigma_n = sigma_n_plus1; 
n = n+1;  
else 
    break
end

end    

sigma_hat = sigma_n;

end

