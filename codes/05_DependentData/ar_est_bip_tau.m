function [phi_hat, x_filt, a_scale_final]=ar_est_bip_tau(x,P)
%   The function ar_est_bip_tau(x,P) comuptes BIP tau-estimates of the
%   AR model parameters. It also computes an outlier cleaned signal using
%   BIP-AR(P) predictions, and the tau-scale of the estimated innovations
%   series.
%
%% INPUTS
% x: data (observations/measurements/signal) 
% P: autoregressive order

%% OUTPUTS
% phi_hat: a vector of bip tau estimates for each order up to P
% x_filt: cleaned version of x using robust BIP predictions
% a_scale_final: minimal tau scale of the innovations of BIP AR or AR
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



%%
N = length(x); % length of the observation vector
kap2 = 0.8724286; %kap=var(eta(randn(10000,1)));
phi_grid = (-.99:0.05:.99); % coarse grid search
fine_grid = (-.99:0.001:.99); % finer grid via polynomial interpolation
a_bip_sc = zeros(size(phi_grid)); % residual scale for BIP-AR on finer grid 
a_sc = zeros(size(phi_grid)); % residual scale for AR on finer grid 


% The following was introduced so as not to predict based on highly contaminated data in the
% first few samples.
x_tran = x(1:min(10,floor(N/2)));
sig_x_tran = madn(x);
x_tran(abs(x_tran)>3*sig_x_tran) = sign(x_tran(abs(x_tran)>3*sig_x_tran))*3*sig_x_tran;
x(1:min(10,floor(N/2))) = x_tran(1:min(10,floor(N/2)));


 if P == 0 % autoregressive order zero
    phi_hat = []; % return empty AR-parameter vector
    a_scale_final(1) = tau_scale(x); % AR(0) residual scale equals observation scale

    
 elseif P == 1 % autoregressive order one
   [x_filt, phi_hat, a_scale_final]=bip_ar1_tau(x,N,phi_grid,fine_grid,kap2); % AR(1) tau-estimates
    
 elseif P > 1 % autoregressive order p   
   [x_filt, phi_hat(1,1), a_scale_final]=bip_ar1_tau(x,N,phi_grid,fine_grid,kap2); % AR(1) tau-estimates
   
    for p=2:P
        for mm=1:length(phi_grid)     
            for pp=1:p-1    
                phi_hat(p,pp) = phi_hat(p-1,pp)-phi_grid(mm)*phi_hat(p-1,p-pp);
            end
            predictor_coeffs = [phi_hat(p,1:p-1)'; phi_grid(mm)]';
            M = length(predictor_coeffs);
            
            if mean(abs(roots([1;predictor_coeffs']))<1)==1
            lambda = ma_infinity(predictor_coeffs, 0, 100); 
            sigma_hat = a_scale_final(1)/sqrt(1+kap2*sum(lambda.^2)); % sigma used for bip-model  
            else 
            sigma_hat = madn(x);
            end
               a=zeros(size(x));
               a2=zeros(size(x));
                for ii=p+1:N   
                    a(ii)=x(ii)-predictor_coeffs*(x(ii-1:-1:ii-M)-a(ii-1:-1:ii-M)+sigma_hat*eta(a(ii-1:-1:ii-M)/sigma_hat));
                    a2(ii)=x(ii)-predictor_coeffs*x(ii-1:-1:ii-M);
                end
            a_bip_sc(mm)=tau_scale(a(p+1:end)); % residual scale for BIP-AR
            a_sc(mm)=tau_scale(a2(p+1:end)); % residual scale for AR     
        end

    poly_approx = polyfit(phi_grid,a_bip_sc,5); % polynomial approximation of residual scale for BIP-AR(p) tau-estimates
    a_interp_scale = (polyval(poly_approx,fine_grid)); % interpolation of  residual scale for BIP-AR(p) tau-estimates to fine grid
    poly_approx2 = polyfit(phi_grid,a_sc,5); % polynomial approximation of  residual scale for AR(p) tau-estimates
    a_interp_scale2 = (polyval(poly_approx2,fine_grid)); % interpolation of  residual scale for AR(p) tau-estimates to fine grid

    [temp ind_max]=min(a_interp_scale);
    phi=-fine_grid(ind_max); % tau-estimate under the BIP-AR(p)
    [temp2 ind_max2]=min(a_interp_scale2);
    phi2=-fine_grid(ind_max2); % tau-estimate under the AR(p)

    % final estimate minimizes the residual scale of the two
     if temp2<temp
        ind_max=ind_max2;
        temp=temp2;
     end

     for pp=1:p-1    
        phi_hat(p,pp) = phi_hat(p-1,pp)-fine_grid(ind_max)*phi_hat(p-1,p-pp);
     end
    phi_hat(p,p) = fine_grid(ind_max);
    
    % final AR(P) tau-scale-estimate depending on phi_hat(p,p)
    if mean(abs(roots([1;phi_hat(p,:)']))<1)==1
    lambda = ma_infinity(phi_hat(p,:), 0, 100); 
    sigma_hat = a_scale_final(1)/sqrt(1+kap2*sum(lambda.^2)); % sigma used for bip-model  
    else 
    sigma_hat = madn(x);
    end
    x_filt = zeros(size(x));
     
        for ii = p+1:N
             a(ii)=x(ii)-phi_hat(p,:)*(x(ii-1:-1:ii-M)-a(ii-1:-1:ii-M)+sigma_hat*eta(a(ii-1:-1:ii-M)/sigma_hat));
        end 
            
     
        for ii = p+1:N
             a2(ii)=x(ii)-phi_hat(p,:)*x(ii-1:-1:ii-M);
        end 
        
        if temp2>temp
            a_scale_final(p+1) = tau_scale(a(p+1:N));
        else
            a_scale_final(p+1) = tau_scale(a2(p+1:N));
        end
        
    end    
    phi_hat = phi_hat(p,:); % BIP-AR(P) tau-estimates
    
        for ii = p+1:N
             x_filt(ii)=x(ii)-a(ii)+sigma_hat*eta(a(ii)/sigma_hat);
        end 
    
 end
 
end

function [x_filt, phi_hat, a_scale_final]=bip_ar1_tau(x,N,phi_grid,fine_grid,kap2)

a_scale_final(1) = tau_scale(x); % AR(0): residual scale equals observation scale
    
     
% grid search for partial autocorrelations   
     for mm = 1:length(phi_grid)
            a = zeros(size(x)); % residuals for BIP-AR
            a2 = zeros(size(x)); % residuals for AR 
            
        lambda = ma_infinity(phi_grid(mm), 0, 100); 
        sigma_hat = a_scale_final(1)/sqrt(1+kap2*sum(lambda.^2)); % sigma used for BIP-model  
         for ii = 2:N
            a(ii) = x(ii)-phi_grid(mm)*(x(ii-1)-a(ii-1)+sigma_hat*eta(a(ii-1)/sigma_hat)); % residuals for BIP-AR
            a2(ii) = x(ii)-phi_grid(mm)*(x(ii-1)); % residuals for AR
         end
    a_bip_sc(mm) = tau_scale(a(2:N)); % tau-scale of residuals for BIP-AR
    a_sc(mm) = tau_scale(a2(2:N)); % tau-scale of residuals for AR
     end

    
    poly_approx = polyfit(phi_grid,a_bip_sc,5); % polynomial approximation of tau scale objective function for BIP-AR(1) tau-estimates
    a_interp_scale = (polyval(poly_approx,fine_grid)); % interpolation of tau scale objective function for BIP-AR(1) tau-estimates to fine grid
    poly_approx2 = polyfit(phi_grid,a_sc,5); % polynomial approximation of  tau scale objective function for AR(1) tau-estimates
    a_interp_scale2 = (polyval(poly_approx2,fine_grid)); % interpolation of  tau scale objective function for AR(1) tau-estimates to fine grid
    
    [temp, ind_max] = min(a_interp_scale);
    phi = fine_grid(ind_max); % tau-estimate under the BIP-AR(1)
    [temp2, ind_max2] = min(a_interp_scale2);
    phi2 = fine_grid(ind_max2); % tau-estimate under the AR(1)

    % final estimate maximizes robust likelihood of the two
     if temp2<temp
        phi_s = phi2;
        temp = temp2;
     else
        phi_s = phi;
     end
    phi_hat = phi_s; % final BIP-tau-estimate for AR(1)
    
    % final AR(1) tau-scale-estimate depending on phi_hat
    lambda = ma_infinity(phi_hat, 0, 100); 
    sigma_hat = a_scale_final(1)/sqrt(1+kap2*sum(lambda.^2));
    a = zeros(size(x)); % residuals for BIP-AR
    a2 = zeros(size(x)); % residuals for AR
    x_filt = zeros(size(x));

        for ii = 2:N
            a(ii) = x(ii)-phi_hat*(x(ii-1)-a(ii-1)+sigma_hat*eta(a(ii-1)/sigma_hat));
            x_filt(ii)=x(ii)-a(ii)+sigma_hat*eta(a(ii)/sigma_hat);
        end 

        for ii = 2:N
            a2(ii) = x(ii)-phi_hat*(x(ii-1));
        end 
        
        if temp2<temp
        a_scale_final(2) = tau_scale(a(2:N));
        else
        a_scale_final(2) = tau_scale(a2(2:N));
        end
        
end