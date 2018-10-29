%%  created by Michael Muma
%   
%   Sets the parameter values for the tracking example.
%
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

    parameter.N   = 3000; % length of simulated trajectory
    parameter.BS  =  [2000 7000; 12000 7000; 7000 12000; 7000 2000; 7000 7000]; % positions of base stations
    parameter.M   = length(parameter.BS(:,1)); % number of base stations
    parameter.start = [4300 4300 2 2]';   % mean starting position
    parameter.pnlos = [ 0 0 0 0 0.25];  
    parameter.initial_sigma = [50 50 6 6]';          % standard deviation for the for the initial state 
    parameter.mc =50 ;        % number of Monte Carlo runs
    parameter.discardN = 100;     % discard first N samples for calculating the error metrics in
    parameter.grid = 1;
    parameter.Ts = 0.2;             %sampling frequency
    parameter.A = [1 0 parameter.Ts 0; 0 1 0 parameter.Ts; 0 0 1 0; 0 0 0 1];
    parameter.G = [parameter.Ts^2/2*eye(2); parameter.Ts*eye(2)];
    parameter.dim = length(parameter.BS(1,:));             % dimension of positions, default is 2 
    parameter.numberbs = 'variable';  
    parameter.noisemodel =  'GMM1'; % noise model
    parameter.motionmodel = 'random-force'; % state model
    parameter.noisestructure = 'Markov'; % noise process model

    
    % graphical output 
    parameter.figure = 1;
    parameter.plot = 'mse';
    
    for ii=1:parameter.M
            switch parameter.pnlos(ii)
                case 0
                    parameter.MarkovChain{ii} = {1 0; 1  0};
                case 0.1
                    parameter.MarkovChain{ii} = {0.99 0.01; 0.09 0.91};
                case 0.25 
                    %parameter.MarkovChain{ii} = {0.98 0.02; 0.06 0.94};%{0.9 0.1; 0.3 0.7};
                    parameter.MarkovChain{ii} = {0.994 0.006; 0.02 0.98};%{0.9 0.1; 0.3 0.7};
                case 0.5
                    %parameter.MarkovChain{ii} = {0.995 0.005; 0.005 0.995};
                    parameter.MarkovChain{ii} = {0.98 0.02; 0.02 0.98};
                case 0.75 
                    parameter.MarkovChain{ii} = {0.94 0.06; 0.02 0.98}; %{0.7  0.3; 0.1 0.9};
                case 1
                    parameter.MarkovChain{ii} = { 0 1; 0 1};
                otherwise
                    disp('Markov Matrix for Probability p_nlos is not available')
            end
    end           
    parameter.sigma_nlos = 400;    % 400, 300 
    parameter.mu_nlos    = 1400;    % 800 1400
    parameter.sigma_v    = 1;
    parameter.sigma_los  = 150;
    parameter.mc=1;
       
    
    
    % parameters for tracker
    % EKF
    ekf.R = parameter.sigma_los^2*diag(ones(1,parameter.M));  
    ekf.Q = parameter.sigma_v*eye(2);
    ekf.G = [parameter.Ts^2/2*eye(2); parameter.Ts*eye(2)];
    ekf.A = parameter.A;
    % initialisation
    ekf.X0 = [30 0 4 10];
    ekf.P0 = diag(parameter.initial_sigma.^2);
    ekf.dim = parameter.dim;
    ekf.var_est = 0;                    
    
    
    
    % Robust EKF 
    rekf.R = parameter.sigma_los^2*diag(ones(1,parameter.M));   
    rekf.Q = parameter.sigma_v*eye(2);
    rekf.G = [parameter.Ts^2/2*eye(2); parameter.Ts*eye(2)];
    rekf.A = parameter.A;
    % initialisation
    rekf.X0 = [30 0 4 10];
    rekf.P0 = diag(parameter.initial_sigma.^2);
    rekf.dim = parameter.dim;
    
    % parameters M-estimator  
    rekf.break    = 1e-4;
    rekf.max_iters = 25;
    rekf.singlescore = 'asymmetric';
    rekf.c1 = 1.5;  % clipping point is c*nu^2
    rekf.c2 = 3.;  % clipping point is c*nu^2
    rekf.var_est = 0;    % uses robust covariance estimation
      x1 =  fzero(@(x1)  rekf.c1-x1*tanh(0.5*x1*(rekf.c2-rekf.c1)),0); %to incorporate nu
    if (x1<0)
        disp('Sign of smoothing parameter is changed')
        rekf.x1 = -x1
    else
        rekf.x1 = x1;
    end
    
    rekf.numberit(parameter.mc,parameter.N) = 0;
    ekf_th(parameter.mc,4,parameter.N) = 0;
    ekf_Hc(parameter.mc,4,parameter.N) = 0;