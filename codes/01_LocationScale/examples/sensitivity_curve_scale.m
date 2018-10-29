
% fix seed of random number generator for reproducibility
rng(2);
N = 100;
x_N_minus1 = randn(N-1,1)+5;
delta_x = linspace(0,10,1000);


% Sensitivity Curve for standard deviation
SC_std = zeros(size(delta_x));
std_hat = std(x_N_minus1);
for ii = 1:length(delta_x)
    SC_std(ii) = N*(std([x_N_minus1; delta_x(ii)])-std_hat);
end


% Sensitivity Curve for median absolute deviation
% that does not coverge to IF
SC_mad = zeros(size(delta_x));
std_hat = madn(x_N_minus1);
for ii = 1:length(delta_x)
    SC_mad(ii) = N*(madn([x_N_minus1; delta_x(ii)])-std_hat);
end

% Sensitivity Curve for mean absolute deviation
% around the median
SC_mead = zeros(size(delta_x));
std_hat = mean(abs(x_N_minus1-median(x_N_minus1)));
for ii = 1:length(delta_x)
    SC_mead(ii) = N*(mean(abs([x_N_minus1; delta_x(ii)]-median(x_N_minus1)))-std_hat);
end

% Sensitivity Curve for Huber's scale estimate
c =  1.3415;
SC_hub = zeros(size(delta_x));
std_hat = MscaleHUB(x_N_minus1,c);
for ii = 1:length(delta_x)
SC_hub(ii) = N*(MscaleHUB([x_N_minus1; delta_x(ii)],c)-std_hat);
end

% Sensitivity Curve for Tukey's scale estimate
c = 4.68; 
SC_tuk = zeros(size(delta_x));
std_hat = MscaleTUK(x_N_minus1,c);
for ii = 1:length(delta_x)
SC_tuk(ii) = N*(MscaleTUK([x_N_minus1; delta_x(ii)],c)-std_hat);
end

figure,
set(gca,'FontSize',18) 
hold on
plot(delta_x,SC_std-min(SC_std),'linewidth',2)
plot(delta_x,SC_mead-min(SC_mead),'linewidth',2)
plot(delta_x,SC_hub-min(SC_hub),'linewidth',2)
plot(delta_x,SC_tuk-min(SC_tuk),'linewidth',2)
grid on
xlabel('Outlier value');
ylabel('Sensitivity curve');
legend('Standard deviation', 'madn', 'Huber M', 'Tukey M');




