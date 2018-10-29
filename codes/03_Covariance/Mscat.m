function [C,invC,iter,flag] = Mscat(X,loss,losspar,invC,printitn)
% [C,invC,iter,flag] = Mscat(X,loss,...)
% computes M-estimator of scatter matrix for the n x p data matrix X  
% using the loss function 'Huber' or 't-loss' and for a given parameter of
% the loss function (i.e., q for Huber's or degrees of freedom v for 
% the t-distribution). 
%
% Data is assumed to be centered (or the symmetry center parameter = 0)
%               
% INPUT:
%        X: the data matrix with n rows (observations) and p columns.
%     loss: either 'Huber' or 't-loss' or 'Tyler'
%  losspar: parameter of the loss function: q in [0,1) for Huber and 
%           d.o.f. v >= 0 for t-loss. For Tyler you do not need to specify
%           this value. Parameter q determines the treshold 
%           c^2 as the qth quantile of chi-squared distribution with p 
%           degrees of freedom distribution (Default q = 0.8). Parameter v 
%           is the def.freedom of t-distribution (Default v = 3)
%           if v = 0, then one computes Tyler's M-estimator
%     invC: initial estimate is the inverse scatter matrix (default = 
%           inverse of the sample covariance matrix) 
% printitn: print iteration number (default = 0, no printing)
% OUTPUT:
%        C: the M-estimate of scatter using Huber's weights
%     invC: the inverse of C
%     iter: nr of iterations
%     flag: flag (true/false) for convergence 
%
% Note: Requires Statistics and Machine Learning Toolbox (calls functions 
% chi2cdf and chi2inv in case of Huber's loss function 
%
% version: Sep 22, 2018
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[n, p] = size(X);

realdata = true;
if ~isreal(X)
    realdata=false; 
end

if nargin < 5
    printitn = 0; 
end

if nargin < 4  || isempty(invC)
    C = X'*X/n; % SCM initial start 
    invC = C\eye(p);
end

switch loss
    case 'Huber'
        ufun = @(t,c) ((t<=c) + (c./t).*(t>c)); % weight function u(t)
        if nargin < 3
            q = 0.90;
        else
            q = losspar;
        end          
        if  isreal(q) && isfinite(q) && (0<q) && (q <1) 
            if realdata
                upar = chi2inv(q,p); % threshold for Huber's weight u(t;.)
                b = chi2cdf(upar,p+2)+(upar/p)*(1-q); % consistency factor
            else              
                upar = chi2inv(q,2*p)/2;  
                b = chi2cdf(2*upar,2*(p+1))+(upar/p)*(1-q); 
            end
        else
            fprintf('input ''losspar''  needs to be a scalar 0<= q < 1\n');
            C = []; invC = []; iter=[];
            return;
        end
        const = (1/(b*n));
    case 't-loss'
        if nargin < 3 
            upar = 3; % d.o.f v=3 is used as the default parameter for t-loss
        else
            upar = losspar; % otherwise use d.o.f. v that was given 
        end        
        if ~isreal(upar) || ~isfinite(upar) || upar < 0
            C = []; invC = []; iter=[];
            error(['Mscat: input ''losspar'' (degrees of freedom ' ...
                'for t-loss needs to be larger than =0 \n']); 
        end  
        if realdata && upar ~=0 
            % this is for real data
            ufun = @(t,v) 1./(v+t); % weight function
            b = tloss_consistency_factor(p,upar);
            const = (upar+p)/(b*n);
        elseif ~realdata && upar ~= 0
           % this is for complex data
           ufun = @(t,v) 1./(v+2*t); % weight function
           b = tloss_consistency_factor(2*p,upar);
           const = (upar+2*p)/(b*n);
        else % the case upar == 0 implies Tyler's M-estimator 
           ufun = @(t,v) 1./t; % weight function
           const = p/n;             
        end      
    case 'Tyler'
           ufun = @(t,v) 1./t; % weight function
           const = p/n;     
    otherwise
        fprintf('please specify loss function as Huber, t-loss or Tyler \n');
        C = []; invC = []; iter=[];
        return;
end

MAX_ITER = 1000; % Max number of iteration
EPS = 1.0e-5;    % Iteration accuracy

iter = 1;
flag = false; % flag for convergence 
while (iter<=MAX_ITER)   
    
      t = real(sum((X*invC).*conj(X),2)); % norms        
      C = const*X'*(X.*repmat(ufun(t,upar),1,p)); 
      d = norm(eye(p)-invC*C,Inf);
      
      if mod(iter,printitn)==0
        fprintf('At iter = %4d, dis=%.6f\n',iter,d);
      end
      
      invC = C\eye(p);
    
      if (d<=EPS) 
         flag = true; 
         break;             
      end         
           
      iter = iter+1;  
end

if (iter==MAX_ITER)
     fprintf(1,'WARNING! Slow convergence: the error of the solution is %f\n',d);  
     fprintf(1,'after %d iterations\n',iter);
end      
end

function b=tloss_consistency_factor(p,v)
% computes the concistency factor b = (1/p) E[|| x ||^2 u_v( ||x||^2)] when
% x ~N_p(0,I). 
sfun = @(x,p,v)  (x.^(p/2)./(v+ x) .* exp(-x/2));
c = 2^(p/2)*gamma(p/2);
q = (1/c)*integral(@(x)sfun(x,p,v),0,Inf);
b = ((v+p)/p)*q; % consistency factor  
end

