function [sig] =madn(y)
% madn computes the normalized median absolute deviation estimate of scale, i.e.,
%
% mu_hat = arg min_mu SUM_i rho_TUK(y_i - mu)
%
%
%   INPUTS: 
%           y: data vector of size N x 1
%
%   OUTPUT:  
%           sig: normalized median absolute deviations scale estimate
% version: Sep 27, 2018
% authors: Michael Muma 

% Compute the auxiliary scale estimate as 
if isreal(y)  
    const =1.4815; 
else
    const =1.20112; 
end 
sig = const*median(abs(y-median(y))); 