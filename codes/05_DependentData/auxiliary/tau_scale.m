function [sigma_hat]= tau_scale(x)
b = 0.398545548533895;  % E(muler_rho2) under the standard normal distribution
sigma_m = m_scale(x);
sigma_hat = sqrt(sigma_m^2/(length(x))*1/b*sum(muler_rho2(x/sigma_m)));
end