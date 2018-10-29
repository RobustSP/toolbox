

% fix seed of random number generator for reproducibility
rng(2);
% number of measurements
N = 100;    
% DC voltage in AWGN
x_N_minus1 = randn(N-1,1)+5;    
% outlier values
delta_x = linspace(0,10,1000);

% sensitivity curve for mean
SC_mean = zeros(size(delta_x));
mu_hat = mean(x_N_minus1);
for ii = 1:length(delta_x)
    SC_mean(ii) = N*(mean([x_N_minus1; delta_x(ii)])-mu_hat);
end

% sensitivity curve for median
SC_med = zeros(size(delta_x));
mu_hat = median(x_N_minus1);
for ii = 1:length(delta_x)
    SC_med(ii) = N*(median([x_N_minus1; delta_x(ii)])-mu_hat);
end

% sensitivity curve for Huber's location estimator
c =  1.3415;
SC_hub = zeros(size(delta_x));
mu_hat = MlocHUB(x_N_minus1,c);
for ii = 1:length(delta_x)
    SC_hub(ii) = N*(MlocHUB([x_N_minus1; delta_x(ii)],c)-mu_hat);
end

% sensitivity curve for Tukey's location estimator
c = 4.68; 
SC_tuk = zeros(size(delta_x));
mu_hat = MlocTUK(x_N_minus1,c);
for ii = 1:length(delta_x)
    SC_tuk(ii) = N*(MlocTUK([x_N_minus1; delta_x(ii)],c)-mu_hat);
end

figure,
set(gca,'FontSize',18) 
hold on
plot(delta_x,SC_mean-mean(SC_mean),'linewidth',2)
plot(delta_x,SC_med-mean(SC_med),'linewidth',2)
plot(delta_x,SC_hub-mean(SC_hub),'linewidth',2)
plot(delta_x,SC_tuk-mean(SC_tuk),'linewidth',2)
grid on
xlabel('Outlier value');
ylabel('Sensitivity curve');
legend('mean', 'median', 'Huber M', 'Tukey M');





