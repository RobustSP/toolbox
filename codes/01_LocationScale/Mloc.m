function [mu_hat] = Mloc(y,lossfun)

% Mloc computes the M-estimates of location using an auxiliary scale
% estimate. It uses the iterative reweighted least squares (IRWLS) algorithm   
%
% INPUTS: 
%       y : (numeric) data vector of size N x 1 
% lossfun : (string) either 'huber' or 'tukey' to identify the desired 
%           loss function
% OUTPUT:  
%           mu_hat: M-estimate of location
% version: Sep 25, 2018
% authors: Michael Muma, Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_hat = Mreg(y,ones(length(y),1),lossfun,median(y));