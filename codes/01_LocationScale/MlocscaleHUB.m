function [mu_hat,sigma_hat] = MlocscaleHUB(y,c)

% Mlocscale computes Huber's M-estimates of location and scale.   
%
%   INPUTS: 
%           y: real valued data vector of size N x 1
%           c: tuning constant c>=0 
%
%   OUTPUT:  
%           mu_hat: Huber's M-estimate of location
%           sigma_hat: Huber's M-estimate of scale
% version: Sep 25, 2018
% authors: Michael Muma, Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 2
   if isreal(y)
       c = 1.3415;
   else
       c = 1.215;
   end % Default: approx 95 efficiency for Gaussian errors
end

[mu_hat, sigma_hat] = hubreg(y,ones(length(y),1),c,madn(y),median(y));