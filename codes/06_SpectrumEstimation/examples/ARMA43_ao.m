% ARMA(4,3) by Moses et al with 2 percent additive outliers

rng(0);

%% Generate Data
N = 512;

% example from Moses
a = [1 -1.316 1.4401 -1.0919 0.83527];

b_0 = 0.13137;
b = [0.13137 0.023543 0.10775 0.03516]/b_0;

p = 4;
q = 3;

w = linspace(0,pi,N/2); 
s = exp(1i*w); 
H = polyval(b,s) ./ polyval(a,s);  

z = randn(N,1);
y = filter(b,a,z);
x = y;  % clean data case

 for jj = 50:51:N
 % 2 percent additive outliers
  y(jj) = y(jj)+1000*randn(1,1);
  end

%% Compute PSD estimates

% periodogram estimate
Pxx = abs(fft(y)).^2;

% BIP Tau ARMA estimate
[PxxdB_tau, Pxx_tau, w, sigma_hat] = spec_arma_est_bip_tau(y,p,q);

% Repeated Median PSD Estimator
[xFRM, ARM, BRM] = repeated_median_filter(y);
Pxx_RM = abs(fft(xFRM)).^2;


%% Plot results
figure, 
set(gca,'FontSize',18) 
hold on
plot(w/(2*pi),10*log10(1/(2*pi)*abs(H).^2), 'linewidth',2) % true spectrum
plot(w/(2*pi), 10*log10(1/N/(2*pi)*Pxx(1:N/2)), 'linewidth',2) % periodogram
xlabel('Normalized frequency')
ylabel('PSD [dB]')
legend('true PSD', 'Periodogram') 


%% plot results
figure, 
set(gca,'FontSize',18) 
hold on
plot(w/(2*pi),10*log10(1/(2*pi)*abs(H).^2), 'linewidth',2) % true spectrum
plot(w/(2*pi),PxxdB, 'linewidth',2) % BIP-ARMA tau PSD estimate
xlabel('Normalized frequency')
ylabel('PSD [dB]')
legend('true PSD', 'BIP tau-estimate') 

%% plot results
figure, 
set(gca,'FontSize',18) 
hold on
plot(w/(2*pi),10*log10(1/(2*pi)*abs(H).^2), 'linewidth',2) % true spectrum
plot(w/(2*pi), 10*log10(1/N/(2*pi)*Pxx_RM(1:N/2)), 'linewidth',2) % repeated median PSD estimate
xlabel('Normalized frequency')
ylabel('PSD [dB]')
legend('true PSD', 'Repeated median estimate') 