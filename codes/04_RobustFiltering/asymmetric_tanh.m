function [phi, phi_point] =  asymmetric_tanh(Sig,c1,c2,x1) 
% The function computes the smoothed assymetric tanh score function and its
% derivative.
%
%% INPUTS
% c1: first clipping point
% c2: second clipping point
% x1 smoothing parameter to make the score function continuous;  %
% x1 =  fzero(@(x1)  c1-x1*tanh(0.5*x1*(c2-c1)),1);

%% OUTPUTS
% phi: assymetric tanh score function
% phi_point: derivative of assymetric tanh score function
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
     
     % phi
    phi(abs(Sig)<= (c1))    = Sig(abs(Sig)<(c1));
    phi(abs(Sig)>c1)        = x1*sign(Sig(abs(Sig)>c1)).*tanh(x1*0.5*(c2-abs(Sig(abs(Sig)>c1))));
    phi(abs(Sig)> c2)       = zeros(1,length(Sig(abs(Sig)>c2)));

    % phi point
    phi_point(abs(Sig)<= (c1)) = 1;
    phi_point(abs(Sig)>c1)  = -0.5*x1^2./(cosh(0.5*x1*(c2-abs(Sig(abs(Sig)>c1))))).^2;
	phi_point(abs(Sig)>c2) = zeros(1,length(Sig(abs(Sig)>c2)));
   
   
return
