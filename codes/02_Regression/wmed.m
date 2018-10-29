function [beta, iter, converged] = wmed(y,w,verbose)
% [beta, iter, converged] = wmed(y,w,...)
% wmed computes the weighted median for data vector y and weights w, i.e. 
% it solves the optimization problem:
%
% beta = arg min_b  SUM_i | y_i - b | * w_i  
%
% inputs: 
%       y : (numeric) data given (real or complex) 
%       w : (nubmer) positive real-valued weights. Inputs need to be of
%           same length
%       verbose: (logical) true of false (default). Set as true if you wish  
%       to see convergence as iterations evolve
% outputs:
%       beta: (numeric) weighted median
%       converged: (logical) flag of convergence (in the complex-valued 
%       data case)
%       iter: (numeric) the number of iterations (complex-valued case)
% version: Aug 31, 2018
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 
    verbose = false;
end
converged = []; iter = [];

N = numel(y); 
if ~isreal(w)
    error('input w needs to a non-negative weight vector'); 
end

if numel(w) ~= N || sum(w >=0) ~= N
    error('wmed: nr of elements of y and w are not equal or w is not non-neg.');   
end

%% real-value case
if isreal(y)
    [y, indx] = sort(y);
    w = w(indx);
    wcum = cumsum(w);
    i = find(wcum<0.5*sum(w),1,'last');
    beta = y(i+1); % due to equation (2.21) of the book  
else 
%% complex-valued case     
    beta0 = median(real(y))+1i*median(imag(y)); % initial guess 
    abs0  = abs(beta0);
    TOL = 1.0e-7;  
    iterMAX = 2000; 
    for iter = 1:iterMAX  
        wy = abs(y- beta0);
        wy(wy <= 10^-6) = 10^-6;
        update = sum(w.*sign(y-beta0))/sum(w./wy);
        beta = beta0 + update;
        delta = abs(update)/abs0;  
        if verbose && mod(iter,10)==0
            fprintf('At iter = %3d, delta=%.8f\n',iter,delta);
        end
        if (delta<=TOL) 
           break;             
        end
        beta0 = beta; abs0 = abs(beta);
    end
    
    if iter == iterMAX
        converged = false;
    else
        converged = true;
    end
    
end