% ANC based on RLS algorithm nomenclature based on slide 23 of 9_Optimum and adaptive filters.pdf
% 
% Usage: [e, new_sig, weight_mat] = RLS(r, noisy_sig, delta, L)
%
% Inputs:
% r  = the vector of desired signal samples of size Ns
% noisy_sig  = the vector of input signal samples of size Ns
% lamda - the weight parameter,
% L = the length (order) of the FIR filter
%
% Outputs:
% e = the output residual error vector of size Ns
% new_sig = output coefficients of noise estimation
% weight_mat = FIR filter parameters

function [err, new_sig, weight_mat] = RLS_method(noisy_sig, R, lamda, L)

    N_sig = length(noisy_sig);
    N_r = length(R);

    % Initialize outputs

    I = eye(L);
    alpha = 0.01;
    p = alpha * I;

    xx = zeros(L,1);                % signal delay matrix
    weight_vect = zeros(L,1);       % adaptive weights vector
    weight_mat = zeros(L,N_sig);    % adaptive weights matrix
    new_sig = zeros(N_sig,1);
    err = zeros(N_sig,1);           % Error vector


    if (N_r <= L)  
        print('error: Signal length is less than the filter order');
        return; 
    end
    if (N_r ~= N_sig)  
        print('error: Input signal and reference signal are different in length?');
        return; 
    end

    for n = 1:N_sig
        xx(1) = R(n);                               % R(n)
        k = (p * xx) ./ (lamda + xx' * p * xx);
        new_sig(n) = xx'*weight_vect;               % new_sig(n) = weight_mat(n)'.R(n)
        err(n) = noisy_sig(n) - new_sig(n);         % e(n) = noisy_sig(n) - weight_mat(n)'.R(n) = noisy_sig(n) - new_sig(n)
        weight_vect = weight_vect + k * err(n);     % weight_mat(n+1) = weight_mat(n) + k**e(n)
        p = (p - k * xx' * p) ./ lamda;
        weight_mat(:,n) = weight_vect;              % Store the weights vector in the weights matrix
        xx(2:L) = xx(1:L-1);                        % Delay back the reference signal window by one sample
    end

end
