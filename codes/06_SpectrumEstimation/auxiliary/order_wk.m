function [wr PSD] = order_wk(y)

% removing a robust location estimate (the sample median) from the
% contaminated time series
ym = median(y);
yt = y - ym;


% initializations for computing an initial robust periodogram, which is
% used to determine the order of fitting the sine and cosine coefficients
% in the repeated median transform

N = length(y);      % length of the time series
K = (N-1)/2;        % number of fourier coefficients
ARM = zeros(1,K);   % repeated median estimate of cosine coefficients at w(k)
BRM = zeros(1,K);   % repeated median estimate of sine coefficients at w(k)
PSD = zeros(1,K);   % robust periodogram for determine the order


for k = 1:K
        w(k) = 2*pi*k/N;
        for u = 0:N-1
            for v = 0:N-1
                if u == v
                    Apuv(u+1,v+1)=0; Bpuv(u+1,v+1)=0;
                else
                    Apuv(u+1,v+1)=( yt(u+1)*sin(w(k)*v)-yt(v+1)*sin(w(k)*u) )./sin(w(k)*(v-u)); % solution to 2 linear equations
                    Bpuv(u+1,v+1)=( yt(v+1)*cos(w(k)*u)-yt(u+1)*cos(w(k)*v) )./sin(w(k)*(v-u));
                end;
            end;
        end;
        Ap(k) = median(median(Apuv')); % repeated median estimate of cosine at frequency w(k)
        Bp(k) = median(median(Bpuv')); % repeated median estimate of sine at frequency w(k)
        PSD(k) = Ap(k).^2+Bp(k).^2;    % robust periodogram to determine the order of fitting the sine and cosine coefficients 
 end;
 
N_win = 3;
Smoothing_Win = hanning(N_win);
PSD_smooth = conv(PSD,Smoothing_Win);
PSD_smooth = PSD_smooth(N_win-1:end-(N_win-2));
 
 % order in which to subtract RM estimates from the time-series is given by
 % I'
 [PSD_sort,I] = sort(PSD',1,'descend');
 order = I';
 wr = 2*pi*order/N;  % fourier frequencies in descending order