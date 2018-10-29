function parameter = create_environment_book(parameter,start,sigma_v)
%%  created by Michael Muma
%   version 30 May 2018
% The function computes generates data for tracking of a mobile user
% equipment based on the values received by parameter. It is based on a function written by Ulrich Hammes, 
% Signal Processing Group, TU Darmstadt, February 2009
%

x = parameter.BS(:,1);
y = parameter.BS(:,2);
parameter.sigma_v = sigma_v;
d(parameter.M,parameter.N)  = 0;
   
            % random-force motion model
            xx(parameter.dim*2,parameter.N) =  0;
            xx(:,1) = start;
            for ii=1:parameter.N-1
                xx(:,ii+1) = parameter.A*xx(:,ii) + parameter.G*parameter.sigma_v *randn(2,1);
            end

            parameter.truetrajectory = xx;
            
            for jj = 1: parameter.M 
                for ii=1: parameter.N
                    d(jj,ii) = sqrt((x(jj)-xx(1,ii))^2+(y(jj)-xx(2,ii))^2);
                end
            end 
            parameter.truedistances = d;
  


        % create noise
        nn(parameter.M,parameter.N)    = 0;
        index(parameter.M,parameter.N) = 0;
   

         for ii=1:parameter.M
            [nn(ii,:), index(ii,:) ] = markov_chain_book(parameter.MarkovChain{ii},parameter.sigma_los,parameter.sigma_nlos,parameter.mu_nlos,parameter.N);
         end

   
            % alignement of variables
            parameter.start             = start;
            size(xx)
            parameter.truetrajectory    = xx;
            parameter.truedistances     = d; 
            parameter.noiseindices      = index;
            parameter.noise = nn;
            parameter.measureddistances = d + nn;


    
    ii = parameter.numbermc;
    parameter.thx(ii,:)  = parameter.truetrajectory(1,:);
    parameter.thy(ii,:)  = parameter.truetrajectory(2,:);
    parameter.thvx(ii,:) = parameter.truetrajectory(3,:);
    parameter.thvy(ii,:) = parameter.truetrajectory(4,:);
    
   
   



return