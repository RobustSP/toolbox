function result = arma_est_bip_mm(x,p,q)
%   The function  arma_est_bip_mm(x,p,q) comuptes BIP MM-estimates of the
%   ARMA model parameters. It also computes an outlier cleaned signal using BIP-ARMA(p,q) predictions
%
%% INPUTS
% x: data (observations/measurements/signal) 
% p: autoregressive order
% q: moving-average order
%
%% OUTPUTS
% result.ar_coeffs: vector of BIP-AR(p) MM-estimates
% result.ma_coeffs: vector of BIP-MA(q) MM-estimates
% result.inno_scale: BIP s-estimate of the innovations scale
% result.cleaned signal: outlier cleaned signal using BIP-ARMA(p,q) predictions
% result.ar_coeffs_init: robust starting point for BIP-AR(p) MM-estimates
% result.ma_coeffs_init: robust starting point for BIP-MA(q) MM-estimates
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



%% Get starting point and residual scale from S-estimtor
bip_s_est = arma_est_bip_s(x,p,q);
beta_hat_s = [bip_s_est.ar_coeffs bip_s_est.ma_coeffs];

%% Solve objective function via nonlinear LS using Levenberg-Marquard algorithm
options = optimset('Display','off','GradObj','on','TolX',5*1e-7,'TolFun',5*1e-7,'LargeScale','on',...
 'Algorithm','Levenberg-Marquardt');

        N = length(x);
        F_mm = @(beta)(1/(N-p)*muler_rho2(arma_s_resid(x, beta, p, q))/bip_s_est.inno_scale);
        F_bip_mm = @(beta)(1/(N-p)*(muler_rho2(bip_s_resid(x, beta, p, q))/bip_s_est.inno_scale));
        
        
        beta_arma_mm = lsqnonlin(F_mm,-beta_hat_s,[],[],options);
        beta_bip_mm = lsqnonlin(F_bip_mm,-beta_hat_s,[],[],options);
        
        a = arma_s_resid(x, beta_arma_mm, p, q);   
        a_sc = m_scale(a); % innovations m-scale for ARMA model
        a_bip = bip_s_resid(x, beta_bip_mm, p, q); 
        a_bip_sc = m_scale(a_bip);% innovations m-scale for BIP-ARMA model
        
        % final parameter estimate uses the model that provides smaller
        % m-scale
        if a_sc<a_bip_sc
            beta_hat = beta_arma_mm;
        else
            beta_hat = beta_bip_mm;
        end
        
        % final m-scale
        a_m_sc = min(a_sc,a_bip_sc);
           
        
%% Output the results
phi_bip_mm = [];
if p>0
    phi_bip_mm = - beta_hat(1:p);
end

theta_bip_mm = [];
if q>0
    theta_bip_mm = -beta_hat(p+1:end);
end

result.ar_coeffs = phi_bip_mm;
result.ma_coeffs = theta_bip_mm;
result.inno_scale = a_m_sc;
result.cleaned_signal = bip_s_est.cleaned_signal;
result.ar_coeffs_init = bip_s_est.ar_coeffs;
result.ma_coeffs_init = bip_s_est.ma_coeffs;   
end        