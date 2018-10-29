function [sigma_hat] = m_scale(x)
N = length(x);
sigma_k =  madn(x);
max_iters = 30;
delta = 3.25/2; % max(muler_rho1)/2
epsilon = 1e-4;
k = 0;
w_k = ones(size(x));

while  and(k<=max_iters,sigma_k<10^5)    
    w_k(x~=0) = muler_rho1(x(x~=0)/sigma_k)./(x(x~=0)/sigma_k).^2;
    w_k(x==0) = 1;
    sigma_k_plus1 = sqrt(1/(N*delta)*sum(w_k.*x.^2)); 

    if abs(sigma_k_plus1/sigma_k-1)>epsilon
    sigma_k = sigma_k_plus1; 
    k = k+1;  
    else 
        break
    end
end    

sigma_hat = sigma_k;
end