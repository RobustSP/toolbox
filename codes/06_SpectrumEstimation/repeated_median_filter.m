function [xFRM, ARM, BRM] = repeated_median_filter(x)
%   The function repeated_median_filter(x) is our implementation of the
%   method described in 
%
%   "High breakdown methods of time series analysis.
%   Tatum, L.G., and Hurvich, C. M.  
%   Journal of the Royal Statistical Society. Series B (Methodological),
%   pp. 881-896, 1993.
%   
%   The code is based on an implementation by Falco Strasser, Signal Processing
%   Group, TU Darmstadt, October 2010.
%
%% inputs
% x: data (observations/measurements/signal), real-valued vector
%
%
%% outputs
% xFRM: repeated median filtered (outlier cleaned) signal 
% ARM:  Fourier coefficients for cosine
% BRM:  Fourier coefficients for sine 
%%  created by Michael Muma
%   version 31 May 2018
%   when using code, cite our work:
%
%   "Robust Statistics for Signal Processing"
%   Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
%   Cambridge University Press, 2018.

N = length(x);    % lenght of signal
xFRM = zeros(size(x));    %  filter cleand signal


% Works for signals of prime length, therefore, signal is split into two
% overlapping segments which are of prime length
x_split = split_into_prime(x);

%% for each prime segment do repeated median filtering
for ii = 1:size(x_split,2)
    
x_part = x_split(:,ii)';

wr = order_wk(x_part);    % Fourier frequencies in descending order

N_prime = length(x_part);      % length of the prime time segment
t = 0:N_prime-1;          % time vector
K = (N_prime-1)/2;        % number of Fourier coefficients

ARM = zeros(1,K);   % repeated median estimate of cosine coefficients at w(k)
BRM = zeros(1,K);   % repeated median estimate of sine coefficients at w(k)

 
% repeated median transform: estimate ARM and BRM starting with stongest
 % w(k), subtract from time series, repeat M times
 
% number of iterations, as recommended in the paper by Tatum and Hurvich
M = 2; 

% remove a robust location estimate (the sample median) 
xm = median(x_part); 
xt = x_part - xm;


for m = 1:M
    for k = 1:K
        for u = 0:N_prime-1
            for v = 0:N_prime-1
                if u == v
                    Auv(u+1,v+1)=0; Buv(u+1,v+1)=0;
                else
                    Auv(u+1,v+1)=( xt(u+1)*sin(wr(k)*v)-xt(v+1)*sin(wr(k)*u) )./sin(wr(k)*(v-u));
                    Buv(u+1,v+1)=( xt(v+1)*cos(wr(k)*u)-xt(u+1)*cos(wr(k)*v) )./sin(wr(k)*(v-u));
                end;
            end;
        end;
        A(k) = median(median(Auv'));
        B(k) = median(median(Buv'));
        xt = xt - A(k)*cos(wr(k)*t) - B(k)*sin(wr(k)*t);
        ARM(k) = ARM(k) + A(k);
        BRM(k) = BRM(k) + B(k);
    end;
end;



%% recover the core process by regression of the repeated median estimates
% onto the independent parameters 
if ii==1
xFRM1 = xm;
 for k = 1:K
    sumAB = ARM(k)*cos(wr(k)*t) + BRM(k)*sin(wr(k)*t);
    xFRM1 = xFRM1 + sumAB;
 end;

elseif ii==2
    xFRM2 = xm;
 for k = 1:K
    sumAB = ARM(k)*cos(wr(k)*t) + BRM(k)*sin(wr(k)*t);
    xFRM2 = xFRM2 + sumAB;
 end
end


end  % end of for ii=1:size(x_split,2) loop

%% fuse the cleaned segments of prime length
if ii==1
   xFRM=xFRM1;
   
elseif ii==2

   xFRM(1:N-N_prime) = xFRM1(1:N-N_prime);
   xFRM(N-N_prime+1:N_prime) = (xFRM1(N-N_prime+1:N_prime) + xFRM2(1:(end-(N-N_prime))))/2;
   xFRM(N_prime+1:end)=xFRM2((end-(N-N_prime)+1):end);
end
   
   
   
   
   
   
   
   
    