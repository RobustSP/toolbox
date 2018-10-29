function [beta_initial, x_filt] = robust_starting_point(x,p,q)
%   The function  robust_starting_point(x,p,q) provides a robust initial estimate for robust ARMA parameter estimation based on BIP-AR(p_long) approximation. It also computes an outlier cleaned signal using BIP-AR(p_long) predictions
%
%% INPUTS
% x: data (observations/measurements/signal) 
% p: autoregressive order
% q: moving-average order
%
%% OUTPUTS
% beta_initial: robust starting point for AR(p)and MA(q) parameters based on BIP-AR(p_long) approximation
% x_filt: outlier cleaned signal using BIP-AR(p_long) predictions
%
%   The function "robust_starting_point" calls "armax" from Matlab System Identification
%   Toolbox to compute classical ARMA parameter estimate based on cleaned
%   data. Replace highlighted code by a different (nonrobust) ARMA parameter estimator if you
%   do not have the toolbox.


if q==0
    p_long = p;
else
    p_long = min(2*(p+q),4); % usually a short AR model provides best results. Change to longer model, if necessary.
end


[~, x_filt, ~] = ar_est_bip_s(x,p_long);


data = iddata(x_filt); % iddata is from System Identification Toolbox 
temp = armax(data,[p q]); % armax is from System Identification Toolbox 
beta_initial = -[temp.a(2:end) temp.c(2:end)]'; 


if or(sum(abs(roots([1 ; -beta_initial(1:p)]))>1),sum(abs(roots([1 ; -beta_initial(p+1:end)]))>1))>0
instable_init = or(sum(abs(roots([1 ; -beta_initial(1:p)]))>1),sum(abs(roots([1 ; -beta_initial(p+1:end)]))>1))>0;
data = iddata(x_filt); % iddata is from System Identification Toolbox 
temp = armax(data,[p q]); % armax is from System Identification Toolbox 
beta_initial = -[temp.a(2:end) temp.c(2:end)]';     
end    

if or(sum(abs(roots([1 ; -beta_initial(1:p)]))>1),sum(abs(roots([1 ; -beta_initial(p+1:end)]))>1))>0
instable_init = or(sum(abs(roots([1 ; -beta_initial(1:p)]))>1),sum(abs(roots([1 ; -beta_initial(p+1:end)]))>1))>0;
[beta_initial, x_filt] = robust_starting_point(x_filt,p,q);
end

if or(sum(abs(roots([1 ; -beta_initial(1:p)]))>1),sum(abs(roots([1 ; -beta_initial(p+1:end)]))>1))>0
instable_init = or(sum(abs(roots([1 ; -beta_initial(1:p)]))>1),sum(abs(roots([1 ; -beta_initial(p+1:end)]))>1))>0;
beta_initial = zeros(size(beta_initial));
end

end