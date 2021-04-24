function [mu_max,R_in,J,e,d_R,d_S,d_E,d_L,d_I,lambda,eta,beta,phi,alpha_S,d_V,d_A,k_A,A0] = return_parameters(input)
    %% parameters assigned
    
    % known values
    mu_max = 1.2;%1.2;
    R_in = 4;%4;
    J = 0;  % 0.2 prev
    e = 5*10^(-7); %5*10^(-7);
    d_R = 0.0;  % 0.2 prev
    d_S = 0.075;
    d_E = 0.075;
    d_L = 0.075;
    d_I = 0.075;
    lambda = 2;
    eta = 1;
    phi = 3.4*10^(-10);
    alpha_S = -0.1;
    
    % doubtful
    d_V = 0.075;              % table says 0.4?
    % unknown
    d_A = 0.1;
    k_A = 7.5*10^7; % 6.022*10^6
    A0 = 10^12; 
    
   % beta = tradeoff_value(input(3));
    beta = 50;
end

