function  [xFB, ABi, BBi] = biweight_filter(x)
%   The biweight_filter(x) is our implementation of the
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
%% INPUTS
% x: data (observations/measurements/signal), real-valued vector
%
%
%% OUTPUTS
% xFBi: Biweight filtered (outlier cleaned) signal 
% ABi:  Fourier coefficients for cosine
% BBi:  Fourier coefficients for sine 
%%  created by Michael Muma
%   version 31 May 2018
%   when using code, cite our work:
%
%   "Robust Statistics for Signal Processing"
%   Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
%   Cambridge University Press, 2018.


N=length(x);    % lenght of signal
xFB=zeros(size(x));    %  filter cleand signal


% Works for signals of prime length, therefore, signal is split into two
% overlapping segments which are of prime length
x_split = split_into_prime(x);


%% for each prime segment do Biweight filtering
for ii=1:size(x_split,2)
    
x_part = x_split(:,ii);
wr = order_wk(x_part);    % Fourier frequencies in descending order
[xFRM, ARM, BRM] = repeated_median_filter(x_part); % initialize with repeated median filter


N_prime = length(x_part);      % length of the prime time segment
t = 0:N_prime-1;          % time vector
K = (N_prime-1)/2;        % number of Fourier coefficients
ABi = zeros(1,K);   % biweight estimate of cosine coefficients at w(k)
BBi = zeros(1,K);   % biweight estimate of sine coefficients at w(k)


k = 4;    % tuning constant as recommended in the paper by Tatum and Hurvich
xb = MlocTUK(x_part, k);   % Tukey's location M-estimate
xc = (x_part-xb)';  % robustly centered time series

for k = 1:K
    s = fitoptions('Method','NonLinearLeastSquares','Robust','Bisquare','Startpoint',[ARM(k);BRM(k)]);
    f = fittype('a*cos(w*x)+b*sin(w*x)','problem','w','options',s) ;
    c = fit(t',xc',f,'problem',wr(k));
    ABi(k) = c.a;
    BBi(k) = c.b;
    xc = xc - ABi(k)*cos(wr(k)*t) - BBi(k)*sin(wr(k)*t);
end


%% recover the core process by regression of the Biweight estimates
% onto the independent parameters 
if ii==1
xFB1 = xb;
    for k = 1:K
        sumAB = ABi(k)*cos(wr(k)*t) + BBi(k)*sin(wr(k)*t);
        xFB1 = xFB1 + sumAB;
    end;

elseif ii==2
    xFB2 = xb;
    for k = 1:K
        sumAB = ABi(k)*cos(wr(k)*t) + BBi(k)*sin(wr(k)*t);
        xFB2 = xFB2 + sumAB;
    end
end


end  % end of for ii=1:size(x_split,2) loop


%% fuse the cleaned segments of prime length
if ii==1
   xFB=xFB1;
   
elseif ii==2

   xFB(1:N-N_prime) = xFB1(1:N-N_prime);
   xFB(N-N_prime+1:N_prime) = (xFB1(N-N_prime+1:N_prime) + xFB2(1:(end-(N-N_prime)))  )/2;
   xFB(N_prime+1:end)=xFB2((end-(N-N_prime)+1):end);
end
   



