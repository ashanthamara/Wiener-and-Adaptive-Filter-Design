function y_hat = weiner_filter(signal, weight_mat)
    order = length(weight_mat);
    y_hat = signal;
    for i = 1: length(signal) - order
        y_hat(i) = signal(i : i + order - 1) * weight_mat;
    end
end