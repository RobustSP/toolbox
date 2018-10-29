function theta_inf = ma_infinity(phi , theta , Q_long) 
theta = theta(:);
phi = phi(:);
Q    =  length(theta);    % MA order 
P    =  length(phi);    % AR order 
theta_inf  =  deconv([1 ; theta ; zeros(Q_long+P+Q,1)] , [1 ; -phi]); 
theta_inf  =  theta_inf(2:Q_long+1); 

if size(theta_inf,2)>size(theta_inf,1) 
   theta_inf  =  theta_inf(:).'; 
end 

end