
function [a_sc] = arma_tau_resid_sc(x, beta_hat, p, q) 

if p>0
phi_hat = beta_hat(1:p);
    if size(phi_hat,2)>size(phi_hat,1)
        phi_hat = phi_hat';
    end
else
    phi_hat = [];
end



if q>0
 theta_hat = beta_hat(p+1:end);
    if size(theta_hat,2)>size(theta_hat,1)
        theta_hat = theta_hat';
    end
else
 theta_hat = [];
end



N = length(x);
r = max(p,q);
a = zeros(size(x));
x_sc = tau_scale(x); 


if r == 0
    
  a_sc = x_sc;  
  a = x;

else   
    
if or(sum(abs(roots([1; -phi_hat]))>1),sum(abs(roots([1; -theta_hat]))>1))
a_sc = 10^10;
else
    
    if p>=1 && q>=1 % ARMA models
        
        for ii=r+1:N 
            
            % ARMA residuals
            a(ii)=x(ii)-phi_hat'*x(ii-1:-1:ii-p)+theta_hat'*a(ii-1:-1:ii-q);

        end 
        
        
    elseif p==0 && q>=1
        
        for ii=r+1:N   % MA models
            
            % MA residuals
            a(ii)=x(ii)+theta_hat'*a(ii-1:-1:ii-q);   
     

        end
        
    elseif p>=1 && q==0 % AR models
        
        for ii=r+1:N
            
            % AR residuals
            a(ii)=x(ii)-phi_hat'*x(ii-1:-1:ii-p);
    

        end
    end


    
	a_sc = tau_scale(a(p+1:end));

end
    
end

end
