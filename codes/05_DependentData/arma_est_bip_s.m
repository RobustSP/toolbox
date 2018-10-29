function result = arma_est_bip_s(x,p,q)
%   The function  arma_est_bip_s(x,p,q) comuptes BIP S-estimates of the
%   ARMA model parameters. It also computes an outlier cleaned signal using BIP-ARMA(p,q) predictions
%
%
%% inputs
% x: data (observations/measurements/signal) 
% p: autoregressive order
% q: moving-average order
%
%% outputs
% result.ar_coeffs: vector of BIP-AR(p) S-estimates
% result.ma_coeffs: vector of BIP-MA(q) S-estimates
% result.inno_scale: BIP s-estimate of the innovations scale
% result.cleaned signal: outlier cleaned signal using BIP-ARMA(p,q) predictions
% result.ar_coeffs_init: robust starting point for BIP-AR(p) S-estimates
% result.ma_coeffs_init: robust starting point for BIP-MA(q) S-estimates
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

%% Making sure that data is in the correct format
if size(x,2) > size(x,1)
x = x';
end

if and(p==0,q==0)
   result.inno_scale = m_scale(x);
   disp('Please choose a nonzero value for p or q') 
elseif length(x)<=(p+q)
   disp('There are too many parameters to estimate for chosen data size. Reduce model order or use a larger data set.') 
   result=struct([]);
else


%% Robust starting point by BIP AR-S approximation
beta_initial = robust_starting_point(x,p,q);

%% Solve objective function via nonlinear LS using Levenberg-Marquard algorithm
options = optimset('Display','off','GradObj','on','TolX',5*1e-7,'TolFun',5*1e-7,'LargeScale','on',...
 'Algorithm','Levenberg-Marquardt');
     
    F = @(beta)(arma_s_resid_sc(x, beta, p, q)); % objective function for ARMA model
    F_bip = @(beta)(bip_s_resid_sc(x, beta, p, q)); % objective function for BIP-ARMA model

    beta_arma = lsqnonlin(F,beta_initial,[],[],options); % Minimize m-scale using ARMA model
    beta_bip = lsqnonlin(F_bip,beta_initial,[],[],options);    % Minimize m-scale using BIP-ARMA model
       
    a_sc = arma_s_resid_sc(x, beta_arma, p, q);   % innovations m-scale for ARMA model
    [a_bip_sc,x_filt] = bip_s_resid_sc(x, beta_bip, p, q); % innovations m-scale for BIP-ARMA model
        
        % final parameter estimate uses the model that provides smaller
        % m-scale
        if a_sc<a_bip_sc
            beta_hat = beta_arma;
        else
            beta_hat = beta_bip;
        end
        
        % final m-scale
        a_m_sc = min(a_sc,a_bip_sc);
        

        
        
%% Output the results
phi_bip_s = [];
phi_bip_s_init = [];
if p>0
    phi_bip_s = - beta_hat(1:p);
    phi_bip_s_init =  - beta_initial(1:p);
end

theta_bip_s = [];
theta_bip_s_init = [];
if q>0
    theta_bip_s = -beta_hat(p+1:end);
    theta_bip_s_init = -beta_initial(p+1:end);
end

result.ar_coeffs = phi_bip_s;
result.ma_coeffs = theta_bip_s;
result.inno_scale = a_m_sc;
result.cleaned_signal = x_filt;
result.ar_coeffs_init = phi_bip_s_init;
result.ma_coeffs_init = theta_bip_s_init;        
end
end