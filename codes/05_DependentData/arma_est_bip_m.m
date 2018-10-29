function result = arma_est_bip_m(x,p,q,beta_hat_s,a_sc_final)
%   The function  arma_est_bip_m(x,p,q) comuptes the BIP M-estimation step for BIP MM estimates of the
%   ARMA model parameters. It can also be used as a stand-alone
%   M-estimator.
%
%% INPUTS
% x: data (observations/measurements/signal) 
% p: autoregressive order
% q: moving-average order
% beta_hat_s: BIP S-estimate
% a_sc_final: M scale estimate of residuals of BIP S-estimate
%
%% OUTPUTS
% result.ar_coeffs: vector of BIP-AR(p) MM-estimates
% result.ma_coeffs: vector of BIP-MA(q) MM-estimates
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



options = optimset('Display','off','TolX',5*1e-5,'TolFun',5*1e-5,'LargeScale','on',...
 'Algorithm','Levenberg-Marquardt');

        N = length(x);
        
        F_mm = @(beta)(sqrt(1/(N-p)*sum(muler_rho2(arma_resid(x/a_sc_final, beta, p, q)))));
        F_bip_mm = @(beta)(sqrt(1/(N-p)*sum(muler_rho2(bip_resid(x/a_sc_final, beta, p, q))/a_sc_final)));
        
        
        beta_arma_mm = lsqnonlin(F_mm,beta_hat_s,[],[],options);
        beta_bip_mm = lsqnonlin(F_bip_mm,beta_hat_s,[],[],options);
        
        a_rho2_mm = 1/(N-p)*sum(muler_rho2(arma_resid(x/a_sc_final, beta_arma_mm, p, q)));
        a_bip_rho2_mm = 1/(N-p)*sum(muler_rho2(bip_resid(x/a_sc_final, beta_bip_mm, p, q)));
        
        if a_rho2_mm<a_bip_rho2_mm
            beta_hat = beta_arma_mm;
        else
            beta_hat = beta_bip_mm;
        end
        
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
end




