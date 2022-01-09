clc;
clear all;

%% load the initial signal
load('idealECG.mat');

ideal_sig = idealECG - mean(idealECG); %making the mean to zero
% sampling frequency 500 Hz
fs = 500;
len = length(ideal_sig);

time  = linspace(0, len-1, len)*(1/fs);

% Data construction (adding noise)
n_wg = awgn(ideal_sig, 10, 'measured'); % adding 10dB guassian noise to the ideal ECG signal

n_50 = 0.2*sin(2*pi*50*time); % creating the other noise signal

noisyECG = n_wg + n_50; % adding the n_50 signal  to the main signal

%% plotting the noise signals
figure('Name','Ideal ECG signal | Fs = 500 Hz'),
plot(time, ideal_sig);
title('Ideal ECG signal | Fs = 500 Hz'), xlabel('Time (s)'), ylabel('Amplitude (mV)');


figure('Name','Noisy ECG signal | Fs = 500 Hz'),
plot(time, noisyECG);
title('Noisy ECG signal | Fs = 500 Hz'), xlabel('Time (s)'), ylabel('Amplitude (mV)');

%% Part 1

% Getting the desired signal

% extracting the ECG pulse from the ideal_sig
time_single_beat = 1.876*fs : 2.043*fs; %84 samples
yi = ideal_sig(time_single_beat);
n = time(time_single_beat);

% extracting the noise from the T wave to P wave
time_noise = 2.032*fs : 2.072*fs; %21 samples
noisy_sig = noisyECG(time_noise);
noise_n = time(time_noise);
noise_est = [noisy_sig noisy_sig noisy_sig noisy_sig]; %make noise estimate same size as desired signal

%% plotting the selected ECG signal and the noise signal

figure('Name','Selected ECG Single Beat and Noise Estimate from T to P wave')

subplot(1,3,1)
plot(n,yi)
xlim([n(1),n(end)])
title('Selected ECG Beat - fs = 500Hz'), xlabel('Time (s)'), ylabel('Amplitude (mV)')

subplot(1,3,2)
plot(linspace(1,21,21),noisy_sig,'r')
title('Noise Estimate - Iso Electric Segment - iso-noise'), xlabel('Time (s)'), ylabel('Amplitude (mV)');

subplot(1,3,3)
plot(linspace(1,84,84),noise_est,'k')
title('Noise Estimate - Iso Electric Segment - Replicated iso-noise'), xlabel('Time (s)'), ylabel('Amplitude (mV)');

%% for an arbitary filter order

order = 15;
noise_added = noisyECG-ideal_sig;
% get the weight matrix
weight_mat = weinerWEIGHTvect(yi, noise_est, order);

%print('Weight matrix of arbitary filer order = %d/n', order);
disp(weight_mat');

% filter the signal with the obtained weight matrix
y_hat = weiner_filter(noisyECG, weight_mat);

figure('Name', 'Weiner Filtering for Arbitary order = 15')
plot(time, ideal_sig, time, noisyECG, time, y_hat)
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG | order = 15')
title('Weiner Filtering noise corrupted ECG M = 15')
xlabel('Time (s)'), ylabel('Voltage (mV)');


%% Find Optimum filter order and the coefficients

ord_range = 55;
mse_mat = NaN(1,ord_range);

for order = 2: ord_range
    weight_mat = weinerWEIGHTvect(yi, noise_est, order);
    y_hat = weiner_filter(noisyECG(time_single_beat),weight_mat); 
    mse_mat(order) = immse(y_hat, ideal_sig(time_single_beat));
end

figure('Name','Weiner Filters MSE vs Different Orders')
plot(mse_mat)
hold on 
[minMSE,minOrd] = min(mse_mat);
scatter(minOrd, minMSE)
title('Weiner Filters MSE vs Different Filter Orders')
xlabel('Filter Order'), ylabel('MSE');

%% Optimum Order
opt_order = minOrd;

opt_weight_mat = weinerWEIGHTvect(yi, noise_est, opt_order);
opt_y_hat = weiner_filter(noisyECG, opt_weight_mat);

figure('Name', 'Weiner Filtering for Optimum order')
plot(time, ideal_sig, time, noisyECG, time, opt_y_hat)
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG | Optimum order')
title('Weiner Filtering noise corrupted ECG with Optimum order')
xlabel('Time (s)'), ylabel('Voltage (mV)');

%% Filter characteristics
fvtool(opt_weight_mat,1)

%% Plotting the spectrum
[Px_noisy_sig, F1_noisy_sig] = periodogram(noisyECG, [], [], fs);
[Px_noise,F2_noise] = periodogram(noise_added, [], [], fs);
[Px_ideal,F3_ideal] = periodogram(ideal_sig, [], [], fs);
[Px_yhat,F4_yhat] = periodogram(opt_y_hat, [], [], fs);

figure('Name','PSD')
semilogy(F1_noisy_sig, Px_noisy_sig, F2_noise, Px_noise, F3_ideal, Px_ideal, F4_yhat, Px_yhat);
legend('Noise corrupted signal','Noise','Desired Signal','Optimum Wiener Filtered Signal')
title('Power Spectral Density'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');

%% Part 2 
% Making a linear model of an ECG

figure('Name','Selected ECG Single Beat');
plot(linspace(1,84,84), yi)
xlim([0,84])
title('Selected ECG Beat - fs = 500Hz'), xlabel('Sample'), ylabel('Amplitude (mV)')


%%
linear_sig = zeros(1,84);
for i = 1: length(yi)

    if i > 79
        linear_sig(i) = 0;
    elseif i > 70 
        linear_sig(i) = ((0.211062+0.108938) - (0))*(i - 79)/(70 - 79) + (0);  
    elseif i > 68
        linear_sig(i) = (0.211062+0.108938);
    elseif i > 55
        linear_sig(i) = ((0.211062+0.108938) - (0))*(i - 55)/(68 - 55) + (0);  
    elseif i > 40
        linear_sig(i) = 0;
    elseif i > 37
        linear_sig(i) = (0 - (-0.548938))*(i - 37)/(40 - 37) + (-0.548938);
    elseif i > 34
        linear_sig(i) = (1.83106 - (-0.548938))*(i - 37)/(34 - 37) + (-0.548938);
    elseif i > 30 
        linear_sig(i) = (1.83106 - (-0.398938))*(i - 30)/(34 - 30) + (-0.398938);
    elseif i > 28
        linear_sig(i) = ((0) - (-0.398938))*(i - 30)/(28 - 30) + (-0.398938);
    elseif i > 16
        linear_sig(i) = 0;
    elseif i > 12
        linear_sig(i) = ((0.0110618+0.128938) - (0))*(i - 16)/(12 - 16) + (0);
    elseif i > 11
        linear_sig(i) = (0.0110618+0.128938);
    elseif i > 8 
        linear_sig(i) = ((0.0110618+0.128938) - (0))*(i - 8)/(11 - 8) + (0);
    else
        linear_sig(i) = 0;
    end
end
linear_sig = linear_sig + ones(1,84)*0.001;
figure('Name','Linear Model | ECG Single Beat');
plot(linspace(1,84,84), linear_sig)
xlim([0,84])
title('Linear ECG Beat'), xlabel('Sample'), ylabel('Amplitude (mV)')

%%
%% for an arbitary filter order

order = 15;
noise_added = noisyECG-ideal_sig;
% get the weight matrix
weight_mat = weinerWEIGHTvect(linear_sig, noise_est, order);

%print('Weight matrix of arbitary filer order = %d/n', order);
disp(weight_mat');

% filter the signal with the obtained weight matrix
y_hat = weiner_filter(noisyECG, weight_mat);

figure('Name', 'Linear Model | Weiner Filtering for Arbitary order = 15')
plot(time, ideal_sig, time, noisyECG, time, y_hat)
legend('Ideal ECG','Noisy ECG','Linear Model | Weiner Filtered ECG | order = 15')
title('Linear Model | Weiner Filtering noise corrupted ECG M = 15')
xlabel('Time (s)'), ylabel('Voltage (mV)');


%% Find Optimum filter order and the coefficients

ord_range = 55;
mse_mat = NaN(1,ord_range);

for order = 2: ord_range
    weight_mat = weinerWEIGHTvect(linear_sig, noise_est, order);
    y_hat = weiner_filter(noisyECG(time_single_beat),weight_mat); 
    mse_mat(order) = immse(y_hat, ideal_sig(time_single_beat));
end

figure('Name','Linear Model | Weiner Filters MSE vs Different Orders')
plot(mse_mat)
hold on 
[minMSE,minOrd] = min(mse_mat);
scatter(minOrd, minMSE)
title('Linear Model | Weiner Filters MSE vs Different Filter Orders')
xlabel('Filter Order'), ylabel('MSE');

%% Optimum Order
opt_order = minOrd;

opt_weight_mat = weinerWEIGHTvect(linear_sig, noise_est, opt_order);
opt_y_hat = weiner_filter(noisyECG, opt_weight_mat);

figure('Name', 'Linear Model | Weiner Filtering for Optimum order')
plot(time, ideal_sig, time, noisyECG, time, opt_y_hat)
legend('Ideal ECG','Noisy ECG','Linear Model | Weiner Filtered ECG | Optimum order')
title('Linear Model | Weiner Filtering noise corrupted ECG with Optimum order')
xlabel('Time (s)'), ylabel('Voltage (mV)');

%% Filter characteristics
fvtool(opt_weight_mat,1)

%% Plotting the spectrum
[Px_noisy_sig, F1_noisy_sig] = periodogram(noisyECG, [], [], fs);
[Px_noise,F2_noise] = periodogram(noise_added, [], [], fs);
[Px_ideal,F3_ideal] = periodogram(ideal_sig, [], [], fs);
[Px_yhat,F4_yhat] = periodogram(opt_y_hat, [], [], fs);

figure('Name','PSD')
semilogy(F1_noisy_sig, Px_noisy_sig, F2_noise, Px_noise, F3_ideal, Px_ideal, F4_yhat, Px_yhat);
legend('Noise corrupted signal','Noise','Desired Signal','Optimum Wiener Filtered Signal')
title('Power Spectral Density'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');


%% 1.2 Freq Derivation

[y_hat_freq, W_f] = weiner_filter_freq(ideal_sig, noise_added, noisyECG);
figure('Name', 'Weiner Filtering Frequency domain')
plot(time, ideal_sig, time, noisyECG, time, y_hat_freq)
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG | Frequency domain')
title('Weiner Filtering Frequency domain'),
xlabel('Time (s)'),ylabel('Voltage (mV)');

%% Plotting the spectrum

[Px_noisy_sig, F1_noisy_sig] = periodogram(noisyECG, [], [], fs);
[Px_noise,F2_noise] = periodogram(noise_added, [], [], fs);
[Px_ideal,F3_ideal] = periodogram(ideal_sig, [], [], fs);
[Px_yhat_freq, F4_yhat_freq] = periodogram(y_hat_freq, [], [], fs);

figure('Name','PSD Weiner Filter Frequency domain')
semilogy(F1_noisy_sig, Px_noisy_sig, F2_noise, Px_noise, F3_ideal, Px_ideal, F4_yhat_freq, Px_yhat_freq);
legend('Noise corrupted signal','Noise','Desired Signal','Frequency domain Wiener Filtered Signal')
title('Power Spectral Density'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');

%% Comparing Frequency domain and Time domain weiner filter

figure('Name', 'Comparing Frequency domain and Time domain weiner filter')
plot(time, ideal_sig, 'g', time, opt_y_hat,'b', time, y_hat_freq, 'r')
legend('Ideal ECG','Filtered by Optimum Time Domain derived Weiner ', 'Filtered by Freq Domain derived Weiner Filter')
title('Comparing Frequency domain and Time domain weiner filter')
xlabel('Time (s)'), ylabel('Voltage (mV)');

mse_time = immse(opt_y_hat, ideal_sig);
mse_freq = immse(y_hat_freq, ideal_sig);

disp('Mean Square Error (Time domain)');
disp(mse_time);
disp('Mean Square Error (Frequency domain)');
disp(mse_freq);


%% 1.3 Effect on non-stationary noise on the Wiener filtering

% creating non stationary noise

f1 = 50;
f2 = 100;

time_p1 = time(1 : floor(length(time)/2));
time_p2 = time(floor(length(time)/ 2) + 1 : end);

n50_p1 = 0.2*sin(2*pi*f1*time_p1);
n100_p1 = 0.3*sin(2*pi*f2*time_p2);

non_stationary_noise = [n50_p1 n100_p1];
non_sta_noisy_sig = n_wg + non_stationary_noise;

% filtering the non stationary noise added signal with the derived frequency domain weiner filter

N = length(non_sta_noisy_sig);  
S_Xf  = fft(non_sta_noisy_sig, N*2-1);
S_Yhat = W_f.* S_Xf;                    % Signal estimate from observation and using Wiener filter
y_hat_time = ifft(S_Yhat);              % converting to time domain
y_hat_non_stat = y_hat_time(1 : N);


figure('Name','Effect of Non Stationary Noise Comparison')
plot(time, ideal_sig, time, y_hat_non_stat, time, y_hat_freq)
xlim([time(1),time(end)])
legend('Ideal ECG','Non-Stationary Noise - Filtered','Stationary Noise - Filtered')
title('Effect of Non Stationary Noise Comparison after filtering with Weiner Freq')
xlabel('Time (s)'), ylabel('voltage(mV)')

%%
noise_added = non_sta_noisy_sig - ideal_sig;

[Px_noisy_sig, F1_noisy_sig] = periodogram(non_sta_noisy_sig, [], [], fs);
[Px_noise,F2_noise] = periodogram(noise_added, [], [], fs);
[Px_ideal,F3_ideal] = periodogram(ideal_sig, [], [], fs);
[Px_yhat_freq_ns, F4_yhat_freq_ns] = periodogram(y_hat_non_stat, [], [], fs);

figure('Name','PSD Weiner Filter Frequency domain | With Non Sattionary Noise')
semilogy(F1_noisy_sig, Px_noisy_sig, F2_noise, Px_noise, F3_ideal, Px_ideal, F4_yhat_freq_ns, Px_yhat_freq_ns);
legend('Non stationary Noise corrupted signal','Non Sattionary Noise','Desired Signal','Frequency domain Wiener Filtered Signal')
title('Power Spectral Density | With Non Sattionary Noise'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');
%%



