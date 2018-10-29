load('ecg_ppg.mat')


%% Plot PPG and ECG (referece) Signals
t_ppg = linspace(0,length(ppg_signal)/1000,length(ppg_signal));
t_ecg = linspace(0,length(ecg_signal)/1000,length(ecg_signal));

figure, 
hold on
subplot(2,1,1); plot(t_ppg,ppg_signal,'linewidth',2);xlabel('Time (s)');ylabel('Amplitude (mV)');title('Photoplethysmogram (PPG)');axis tight;
subplot(2,1,2); plot(t_ecg,ecg_signal,'linewidth',2);xlabel('Time (s)');ylabel('Amplitude (mV)');title('Electrocardiogram (ECG)');axis tight;


%% Plot inter-beat-intervals
ibi_ppg = diff(ppg_pos);
ibi_ecg = diff(ecg_pos);

figure, 
subplot(2,1,1); plot(ibi_ppg,'linewidth',2);xlabel('Samples');ylabel('Inter-beat-interval (ms)');title('Inter-beat-intervals for Photoplethysmogram (PPG)');axis tight;
subplot(2,1,2); plot(ibi_ecg,'linewidth',2);xlabel('Samples');ylabel('Inter-beat-interval (ms)');title('Inter-beat-intervals for Electrocardiogram (ECG)');axis tight;



%% BIP-Tau Data Cleaning of PPG inter-beat-intervals


% Order Selection was computed offline based on robust information criteria
% for pp = 0:10
%     for qq = 0:10
%         result = arma_est_bip_tau(ibi_ppg-median(ibi_ppg),pp,qq);
%         SIC(pp+1,qq+1) = log(result.inno_scale^2)+(pp+qq)/length(ibi_ppg)*log(length(ibi_ppg));
%         AIC(pp+1,qq+1) = log(result.inno_scale^2)+2*(pp+qq)/length(ibi_ppg);
%         HQC(pp+1,qq+1) = log(result.inno_scale^2)+2*(pp+qq)/length(ibi_ppg)*log(log(length(ibi_ppg)));
%     end
% end

% Selected ARMA model
p = 0;
q = 11;

% BIP-ARMA parameter estimation and data cleaning 
result = arma_est_bip_tau(fliplr(ibi_ppg)-median(ibi_ppg),p,q); % data is flipped since outliers appear in the first few samples
ibi_ppg_cl = fliplr(result.cleaned_signal') + median(ibi_ppg);

disp('Selected ARMA model: MA(11)')
disp('estimated coefficients:')
ma_coeffs = result.ma_coeffs
figure, 
subplot(2,1,1); plot(ibi_ppg); hold on; plot(ibi_ppg_cl,'linewidth',2);xlabel('Samples');ylabel('Inter-beat-interval (ms)');title('Cleaned Inter-beat-intervals for Photoplethysmogram (PPG)');legend('original','cleaned');axis tight;
subplot(2,1,2); plot(ibi_ecg); hold on; plot(ibi_ppg_cl,'linewidth',2);xlabel('Samples');ylabel('Inter-beat-interval (ms)');title('Cleaned Inter-beat-intervals for PPG and ECG reference');legend('ECG reference','cleaned');axis tight;
