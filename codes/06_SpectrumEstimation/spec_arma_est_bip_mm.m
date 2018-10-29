function [PxxdB, Pxx, w, sigma_hat] = spec_arma_est_bip_mm(x,p,q)
%   The function spec_arma_est_bip_tau(x,p,q) comuptes spectral estimates using the BIP MM-estimates of the
%   ARMA model parameters.
%
%% INPUTS
% x: data (observations/measurements/signal) 
% p: autoregressive order
% q: moving-average order
%
%% OUTPUTS
% PxxdB: spectral estimate in dB
% Pxx: spectral estimate
% w: frequency (0,pi)
% sigma_hat: BIP M-scale estimate of the innovations
%
%   The function "robust_starting_point" calls "armax" from Matlab System Identification
%   Toolbox to compute classical ARMA parameter estimate based on cleaned
%   data. Replace highlighted code by a different (nonrobust) ARMA parameter estimator if you
%   do not have the toolbox.
%
%   created by Michael Muma
%   version 02 October 2018
%   when using code, cite our work:
%
%   "Robust Statistics for Signal Processing"
%   Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
%   Cambridge University Press, 2018.
%
%   and 
%
%  "Bounded Influence Propagation $\tau$-Estimation: A New Robust Method for ARMA Model Estimation." 
%   Muma, M. and Zoubir, A.M.
%   IEEE Transactions on Signal Processing, 65(7), 1712-1727, 2017.
%


x = x - median(x);
N = length(x);
w = linspace(0,pi,N/2); 
s = exp(1i*w); % Digital frequency must be used for this calculation

result = arma_est_bip_mm(x,p,q);
beta_hat = [result.ar_coeffs; result.ma_coeffs];
Xx = polyval([1; beta_hat(p+1:end)],s) ./ polyval([1; beta_hat(1:p)],s); 
sigma_hat = result.inno_scale;
Pxx = sigma_hat^2/(2*pi)*abs(Xx).^2; % BIP-ARMA tau spectral estimate
PxxdB = 10*log10(Pxx); % BIP-S ARMA spectral estimate in dB

end