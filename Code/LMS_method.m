% ANC based on LMS algorithm nomenclature based on slide 15 of 9_Optimum and adaptive filters.pdf
% 
% Usage: [err_vect, new_sig, weight_mat] = LMS(Rr, noisy_sig, delta, L)
%
% Inputs:
% Rr  = the vector of desired signal samples of size Ns
% noisy_sig  = the vector of input signal samples of size Ns
% miu = rate of convergence
% L = the length (order) of the FIR filter
%
% Outputs:
% err_vect = the output residual error vector of size Ns
% new_sig = output coefficients of noise estimation
% weight_mat = FIR filter parameters

function [err_vect, new_sig, weight_mat] = LMS_method(noisy_sig, sig_R, mu, M)

    N_r = length(sig_R);
    N_x = length(noisy_sig);
    
    % Initialize outputs and required matrices
    Rr = zeros(M,1);                % signal delay matrix
    weight_vect = zeros(M,1);       % adaptive weights vector
    weight_mat = zeros(N_x,M);      % adaptive weights matrix
    new_sig = zeros(1, N_x);
    err_vect = zeros(1, N_x);       % Error vector

    if (N_r <= M)  
        disp('error: Signal length is less than the filter order');
        return; 
    end
    if (N_r ~= N_x)  
        disp('error: Input signal and reference signal are different in length?');
        return; 
    end

    lambda_max = 20*M*((noisy_sig*noisy_sig')/length(noisy_sig));
    
    if (mu > 2/lambda_max) 
        disp(['mu is too large' num2str(mu) ' /' num2str(lambda_max)]);
        return
    end


    for k = 1:N_x
        Rr(1) = sig_R(k);                                   % sig_R(n)
        new_sig(k) = weight_vect'*Rr;                       % new_sig(n) = weight_mat(n)'.sig_R(n) new_sig gives the error
        err_vect(k) = noisy_sig(k) - new_sig(k);            % err_vect(n) = noisy_sig(n) - weight_mat(n)'.sig_R(n) = noisy_sig(n) - new_sig(n)
        weight_vect = weight_vect + 2*mu*err_vect(k)*Rr;    % weight_mat(n+1) = weight_mat(n) + 2*miu*err_vect(n)*sig_R(n) - Widrow Hoff LMS algorthm
        weight_mat(k,:) = weight_vect;                      % Store the weights vector in the weights matrix
        Rr(2:M) = Rr(1:M-1);                                % Delay back the reference signal window by one sample
    end

end