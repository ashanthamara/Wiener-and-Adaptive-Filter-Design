clc;
clear all;

%% Data Construction
fs = 500;

N = 5000;                                                   % Number of samples
time = linspace(0, 5, N);                                   % Time vector
sawtooth_sig = sawtooth(2*pi*2*time(1 : N), 0.5);           % Sawtooth signal
n50 = 0.2*sin(2*pi*50*time(1 : N/2));                       % Sinusoidal noise with 50 Hz
n100 = 0.3*sin(2* pi*100*time(N/2 + 1 : N));                % Sinusoidal noise with 100 Hz
nwg = sawtooth_sig - awgn(sawtooth_sig, 10, 'measured');    % 10 dB Gaussian white noise

noisy_signal = sawtooth_sig + nwg + [n50 n100];             % corrupted signal with noise 

figure('Name', 'Signals using for Adaptive Filtering')
subplot(2,1,1)
plot(time, sawtooth_sig)
title('Sawtooth wave (ideal signal)')
subplot(2,1,2)
plot(time, noisy_signal)
title('Noise added signal (noisy signal)')
linkaxes()

%% Create the signal sig_R

% arbitary constants
a = 1.611; 
phi_1 = pi*(1/6); 
phi_2 = pi*(1/2);

n50_r = 0.2*sin(2*pi*50*time(1 : N/2) + phi_1); 
n100_r = 0.3*sin(2* pi*100*time(N/2 + 1 : N) + phi_2); 

sig_R = a*(nwg + [n50_r n100_r]);

%% LMS_method Algorithm
%% 2.1 a

%%
mu1 = 0.006112;
M = 12;
[err, ~, ~] = LMS_method(noisy_signal, sig_R, mu1, M);
figure;
subplot(4,1,1)
plot(time,sawtooth_sig)
title('Desired Signal sawtooth_sig(n)')
subplot(4,1,2)
plot(time,noisy_signal)
title('Noise Corrupted Signal')
subplot(4,1,3)
plot(time,err)
title(['Filtered Signal using LMS_method M=10 u=0.006161'])
subplot(4,1,4)
plot(time,abs(err-sawtooth_sig))
title('Absolute Error')

e_LMS = err;

%%
MRange = 15;
mse = NaN(MRange, 100);

lambda_max = 20*MRange*((noisy_signal*noisy_signal')/ length(noisy_signal));
mu = linspace(0.001, 2/ lambda_max, 100);

for M = 1:MRange
    for i = 1:100
        [err, ~, ~] = LMS_method(noisy_signal, sig_R, mu(i), M);
        mse(M,i) = immse(err,sawtooth_sig);
    end
end

%%
M = 1:MRange;
surf(mu, M, mse)
title('MSE Variation when Adaptive Filtering using LMS') 
colorbar
xlabel('mu'), ylabel('M - Order'),zlabel('MSE');
colormap('jet')

%% Find LMS with min MSE
[ms, ls] = min(mse,[],2);
[mse_min, m_min] = min(ms);
lambda_min = ls(m_min)*(2/ lambda_max - 0.001)/100 + 0.001;
disp(['Minimum Error = ' num2str(mse_min) ' at M = ' num2str(m_min) ' and mu = ' num2str(lambda_min)])

%%  2.2 RLS Algorithm

%% fastest RLS_method filter
lamda = 0.996;          % 0 < lamda <= 1
M = 15;
[err2, ~, ~] = RLS_method(noisy_signal, sig_R, lamda, M);

figure;
subplot(4,1,1)
plot(time, sawtooth_sig)
title('Desired Signal sawtooth_sig(n)')
subplot(4,1,2)
plot(time, noisy_signal)
title('Noise Corrupted Signal')
subplot(4,1,3)
plot(time, err2);
title('Filtered Signal using RLS_method M=15 lambda = 0.996');
subplot(4,1,4)
plot(time, abs(err2 - sawtooth_sig'))

e_RLS = err2;

%%  
MRange = 15;
mse = NaN(MRange,100);
lambda = linspace(0.9,1,100);

for M=1:MRange
    for i = 1:100
        err = RLS_method(noisy_signal, sig_R, lambda(i), M);
        mse(M,i) = immse(err', sawtooth_sig);
    end
end

%%
figure;
surf(lambda,(1:MRange), mse)
colorbar
title('MSE Variation when Adaptive Filtering using RLS')
xlabel('lambda'), ylabel('M - Order'), zlabel('MSE');
colormap('jet')

%% Find RLS with min MSE
[ms,ls] = min(mse,[],2);
[mse_min,m_min] = min(ms);
lambda_min = ls(m_min)*(0.01)/100 + 0.9;
disp(['Minimum Error = ' num2str(mse_min) ' at M = ' num2str(m_min) ' and lambda = ' num2str(lambda_min)])

%% Comparing LMS_method and RLS_method

figure
subplot(2,1,1)
plot(time, abs(e_LMS - sawtooth_sig));
title(['Error convergence using the LMS algorithm \mu = ' num2str(mu1) ' M = 15' ]);
xlabel('Time(sawtooth_sig)');
ylabel('Voltage (mV)');
grid on
subplot(2,1,2)
plot(time, abs(e_RLS' - sawtooth_sig))
title(['Filtered Signal of the ANC filter using the RLS algorithm \lambda = ' num2str(lamda) ' M = 15']),
xlabel('Time(sawtooth_sig)');
ylabel('Voltage (mV)');
grid on


%% c. Adaptive FIltering ECG_sig signal
load('idealECG.mat')
ECG_sig = idealECG - mean(idealECG);

fs = 500;
N = length(ECG_sig);                                    % Number of points
time = linspace(0,N/fs,N);                  

signal = ECG_sig;                                 
n50 = 0.2*sin(2*pi*50*time(1 : N/2));                     
n100 = 0.3*sin(2* pi*100*time(N/2+1 : N));                

nwg = signal - awgn(signal, 10,'measured');  % Gaussian white noise

noisy_signal = signal + nwg + [n50 n100];         % noise corrupted signal


%% Generating sig_R
a = 1.611; 
phi_1 = pi*(1/6); 
phi_2 = pi*(1/2);

n50_r = 0.2*sin(2*pi*50*time(1 : N/2) + phi_1); 
n100_r = 0.3*sin(2*pi*100*time(N/2 + 1 : N) + phi_2); 

sig_R = a*(nwg + [n50_r n100_r]);

%% LMS to ECG_sig
mu = 0.006112;
M = 12;
[err, ~, ~] = LMS_method(noisy_signal, sig_R, mu, M);

%% RLS to ECG_sig
lamda = 0.996; % 0 < lamda <= 1
M = 12;
[err2, y2, w2] = RLS_method(noisy_signal, sig_R, lamda, M);

%%
figure;
subplot(6,1,1)
plot(time, signal)
title('Desired Signal sawtooth signal(n)')
subplot(6,1,2)
plot(time, noisy_signal)
title('Noise Corrupted Signal')
subplot(6,1,3)
plot(time, err)
title(['Filtered Signal using LMS method M=12 u=0.006112'])
subplot(6,1,4)
plot(time, abs(err - signal))
title('Absolute Error | LMS')
subplot(6,1,5)
plot(time,err2);
title(['Filtered Signal using RLS method M=12 lambda = 0.996']);
subplot(6,1,6)
plot(time,abs(err2 - signal'))
title('Absolute Error | RLS')
linkaxes()