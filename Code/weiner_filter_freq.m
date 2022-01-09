function [y_hat, W_f] = weiner_filter_freq(template, noise, signal)

    N = length(signal);               % All FFT all has to be same size
    
    S_yy = abs(fft(template,N*2-1)).^2;  % Power of template
    S_NN = abs(fft(noise,N*2-1)).^2;     % Power of noise
    S_Xf  = fft(signal,N*2-1);           % FT of the signal
    
    W_f = S_yy./(S_yy + S_NN);          % Freq Wiener filter 
    S_Yhat = W_f.*S_Xf;                 % Signal estimate from observation and using Wiener filter
    
    y_hat_time = (ifft(S_Yhat));            % time domain
    y_hat = y_hat_time(1: length(signal));   % Since FFT being two sided

end