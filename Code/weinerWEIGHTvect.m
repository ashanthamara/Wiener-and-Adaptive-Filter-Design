function W = weinerWEIGHTvect(y, n_est, order)
    % For autocorrelation
    yyT = 0; 
    nnT = 0;
    % For crosscorrelation
    Yy = 0;

    y_mat = zeros(order,1);
    n_mat = zeros(order,1);
    
    len = length(y);

    for i=1:len
        
        y_mat(1) = y(i); 
        n_mat(1) = n_est(i);

        yyT = yyT + toeplitz(autocorr(y_mat, order-1));
        nnT = nnT + toeplitz(autocorr(n_mat, order-1));
        %disp(yyT);

        Yy = Yy + y_mat*y(i);

        % shifting the delay 
        y_mat(2:order) = y_mat(1 : order-1);
        n_mat(2:order) = n_mat(1 : order-1);
    end

    yyT = yyT.*mean(y.^2);
    nnT = nnT.*mean(n_est.^2);

    autocorr_Y = yyT./ (len - order);
    autocorr_N = nnT./ (len - order);
    theta_Yy = Yy./ (len-order);
    
    autocorr_X = autocorr_Y + autocorr_N;
    W = autocorr_X\theta_Yy;
    %W = inv(autocorr_X)*theta_Yy;

end