%%  created by Michael Muma
%   version 30 May 2018
%   when using code, cite our work:
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
%
%
%   The example uses functions that were written by Ulrich Hammes, Signal Processing Group, TU Darmstadt, February 2009
%
%

clear all;
 
% set the parameters for data generation and for trackers
set_parameters_book
        
tic   
for ii=1:parameter.mc 
      
    parameter.numbermc = ii;
    
    % generate measurements
    parameter = create_environment_book(parameter,parameter.start,parameter.sigma_v);
     
    % generate random starting point
    randnvector = [ parameter.initial_sigma(1)*randn parameter.initial_sigma(2)*randn parameter.initial_sigma(3)*randn parameter.initial_sigma(4)*randn]';  
    theta_init = parameter.start + randnvector;
    
    % estimate positions using (robust) extended Kalman filter
    [ekf_th(ii,:,:)]   = ekf_toa(parameter.measureddistances,theta_init,parameter.BS,ekf);
    [ekf_Hc(ii,:,:)]   = ekf_toa_robust(parameter.measureddistances,theta_init,parameter.BS,rekf);
end 
      
toc

show_results_book