% Examples for Bounded Innovation Propagation (BIP) Tau ARMA parameter estimation
rng(0);


%% Example 1: AR(1) with 30 percent isolated outliers

% Generate AR(1) observations
N = 300; a = randn(N,1); x = filter(1,[1 -0.8],a); p = 1; q = 0;

% Generate isolated outliers
cont_prob = 0.3;    % outler contamination probability
outlier_ind = find(sign(rand(N,1)-(cont_prob))<0); % outlier index
outlier = 100*randn(N,1); % contaminating process
v = zeros(N,1); % additive outlier signal
v(outlier_ind) = outlier(outlier_ind);
v(1) = 0; % first sample should not be an outlier

% 30 percent of isolated additive outliers
x_ao = x+v;

% BIP Tau estimation
result = arma_est_bip_tau(x_ao,p,q);
disp('Example: AR(1) with ar_coeff = -0.8')
disp('30 percent isolated additive outliers')
disp('estimated coefficients:')
ar_coeff_est = result.ar_coeffs

figure, 
hold on
subplot(2,1,1);  plot(x_ao); hold on; plot(result.cleaned_signal,'linewidth',2);xlabel('Samples');ylabel('Amplitude');title('BIP-AR(1) cleaned signal');legend('outlier contaminated AR(1)','cleaned');axis tight;
subplot(2,1,2);  plot(x); hold on; plot(result.cleaned_signal,'linewidth',2);xlabel('Samples');ylabel('Amplitude');title('BIP-AR(1) cleaned signal');legend('original AR(1)','cleaned');axis tight;




%% Example 2: ARMA(1,1) with 10 percent patchy outliers

% Generate ARMA(1,1) observations
N = 1000; a = randn(N,1); x = filter([1 0.2],[1 -.8],a); p = 1; q = 1;

% Generate a patch of outliers of length 101 samples
v = 1000*randn(101,1);

% 10 percent of patchy additive outliers
x_ao = x; x_ao(100:200) = x(100:200)+v; 

% BIP Tau estimation
result = arma_est_bip_tau(x_ao,p,q);

disp('Example 2: ARMA(1,1) with ar_coeff = -0.8, ma_coeff 0.2')
disp('10 percent patchy additive outliers')
disp('estimated coefficients:')
ar_coeff_est = result.ar_coeffs
ma_coeff_est = result.ma_coeffs

figure, 
hold on
subplot(2,1,1);  plot(x_ao); hold on; plot(result.cleaned_signal,'linewidth',2);xlabel('Samples');ylabel('Amplitude');title('BIP-ARMA(1,1) cleaned signal');legend('outlier contaminated ARMA(1,1)','cleaned');axis tight;
subplot(2,1,2);  plot(x); hold on; plot(result.cleaned_signal,'linewidth',2);xlabel('Samples');ylabel('Amplitude');title('BIP-ARMA(1,1) cleaned signal');legend('original ARMA(1,1)','cleaned');axis tight;





%% Example 3: MA(2) with 20 percent isolated outliers

% Generate MA(2) observations
N = 500; a = randn(N,1); x = filter([1 -0.7 0.5],1,a); p = 0; q = 2;

% Generate isolated outliers
cont_prob = 0.2;
outlier_ind = find(sign(rand(N,1)-(cont_prob))<0);
outlier = 100*randn(N,1);
v = zeros(N,1);
v(outlier_ind) = outlier(outlier_ind);
v(1:2) = 0;

% 20 percent of isolated additive outliers
x_ao = x+v;

% BIP Tau estimation
result = arma_est_bip_tau(x_ao,p,q);

disp('Example 3: MA(2) with ma_coeff = [-0.7 0.5]')
disp('10 percent patchy additive outliers')
disp('estimated coefficients:')
ma_coeff_est = result.ma_coeffs


figure, 
hold on
subplot(2,1,1);  plot(x_ao); hold on; plot(result.cleaned_signal,'linewidth',2);xlabel('Samples');ylabel('Amplitude');title('BIP-MA(2) cleaned signal');legend('outlier contaminated MA(2)','cleaned');axis tight;
subplot(2,1,2);  plot(x); hold on; plot(result.cleaned_signal,'linewidth',2);xlabel('Samples');ylabel('Amplitude');title('BIP-MA(2) cleaned signal');legend('original MA(2)','cleaned');axis tight;






