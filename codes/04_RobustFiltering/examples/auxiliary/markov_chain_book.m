function [y, h] = markov_chain_book(MC,sigma_los,sigma_nlos,mu_nlos,N)

  MM = cell2mat(MC);
  % transition probabilities 
  pt_losnlos = MM(1,2);     % from LOS->NLOS 
  pt_nloslos = MM(2,1);     % from NLOS->LOS 
  
    m=1;
    h(N+100) = 0;
    y(N+100) = 0;
    pp =rand;

%initialisation
if (pp<0.5)
   
    state=0;    % LOS
    
else
    
    state=1;    % NLOS
end


while m<=N+100
    
    switch state
    
        case 0      % LOS
            
            p=rand(1);
		 
          if (m>100)        % ensures that the Markov Chain is in a steady state (p_(k) = p_init*A^k), if k large, p(k) is steady, where A is the Markov transition matrix  
                y(m) = sigma_los*randn;
                h(m) = 0;
          end
            if p<pt_losnlos
                
                state=1;
                
            end
            
		case 1      % NLOS
            
             q=rand(1);
			 
          if (m>100)	  % ensures that the Markov Chain is in a steady state (p_(k) = p_init*A^k), if k large, p(k) is steady, where A is the Markov transition matrix 
                y(m) = sigma_los*randn + mu_nlos + sigma_nlos*randn;
                h(m) = 1;
          end
             if q<pt_nloslos
                 
                 state=0;
                 
             end
             
    end
    
m=m+1;

end

% the first 100 samples are discarded to ensure that the Markov Chain is in
% the steady state
y = y(101:end);
h = h(101:end);

return