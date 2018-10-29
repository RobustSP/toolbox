function [th2, Theta, kk, residuals] = m_param_est(receive,C,Theta,parameter)
% The function computes an M-estimator of regression using the assymetric tanh score function. 
%
%% INPUTS
% recieve: is the received signal
% Theta: contains the initial estimate
% C: is the regression matrix of the model y = C*x + n 
% parameter.maxiters: maximal number of iterations
% parameter.break: break condition 

%% OUTPUTS
% th2: M-estimate of regression
% Theta:  M-estimates of regression for all iterations
% kk: iteration index
% residual: residuals given Theta
%%  created by Michael Muma
%   Based on a function written by Ulrich Hammes, Signal Processing Group, TU Darmstadt, February 2009
%   version 30 May 2018
%   When using code, cite our work:
%
%   "Robust Statistics for Signal Processing"
%   Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
%   Cambridge University Press, 2018.
%
%   and 
%
%  "Robust Tracking and Geolocation for Wireless Networks in NLOS Environments." 
%   Hammes, U., Wolsztynski, E., and Zoubir, A.M.
%   IEEE Journal on Selected Topics in Signal Processing, 3(5), 889-901, 2009.

    Pseudoinv = pinv(C);
    residuals = receive - C*Theta;
    Theta     = [zeros(length(Theta(:,1)),1)  Theta zeros(length(Theta(:,1)),parameter.max_iters-1)];
    kk = 1;
    
    noisescale = 1.483*mad(residuals);
    
    

 while  (sum(abs((Theta(:,kk+1)-Theta(:,kk))./Theta(:,kk+1)))>parameter.break) % (sum(abs((Theta(:,kk+1)-Theta(:,kk))))>parameter.break)
      
         [z, Phi_point]  =  asymmetric_tanh(residuals'/noisescale,parameter.c1,parameter.c2,parameter.x1);


        z  =   z';
        
        kk = kk + 1;
         
        muu = 1.25* max(abs(Phi_point));
        
        % update parameter estimate
        Theta(:,kk+1) = Theta(:,kk) + 1/muu * Pseudoinv *z*noisescale;
       
        
        residuals = receive - C*Theta(:,kk+1);
        noisescale = 1.483*mad(residuals);
    
        if (kk>=parameter.max_iters)
             break;
         end
  
 end

th2 = Theta(:,kk+1);



return

