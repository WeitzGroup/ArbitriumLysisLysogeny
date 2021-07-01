function dfdtheta = compute_dfdtheta(x_hat, theta, noise_on, fixed_prob, pars)
    
    %% Assign letters from x_hat
    R = x_hat(1);
    S = x_hat(2);
    E = x_hat(3);
    L = x_hat(4);
    I = x_hat(5);
    V = x_hat(6);
    A = x_hat(7);

    
    %% compute monod and probability
    Psi = monod(R,pars.mu_max,pars.R_in);
    P_A = probability(A,theta,pars.A0, noise_on, fixed_prob);
    dfdtheta = zeros(7,2);
    
    %% compute needed derivatives
    
    dPsidR = pars.mu_max/(R+pars.R_in) - pars.mu_max*R/((R+pars.R_in)^2);
    
    max_val = 1; k = theta(1); x0 = theta(2);
    logistic_fn = 1/(1+exp(-k*(A/(pars.A0*x0)-1)));
    dP_Adk  = max_val*(A/(pars.A0*x0)-1)*(logistic_fn*(1-logistic_fn));
    dP_Adx0 = -max_val*k*A/(pars.A0*x0^2)*(logistic_fn*(1-logistic_fn)); 
      
    %% Assign value to dfdx term by term
    
    dfdtheta(1,1) = 0;
    dfdtheta(1,2) = 0;
    
    dfdtheta(2,1) = 0;
    dfdtheta(2,2) = 0;
    
    dfdtheta(3,1) = 0;
    dfdtheta(3,2) = 0;
    
    dfdtheta(4,1) = pars.lambda*E*dP_Adk;
    dfdtheta(4,2) = pars.lambda*E*dP_Adx0;
    
    dfdtheta(5,1) = -pars.lambda*E*dP_Adk;
    dfdtheta(5,2) = -pars.lambda*E*dP_Adx0;
    
    dfdtheta(6,1) = 0;
    dfdtheta(6,2) = 0;
    
    dfdtheta(7,1) = 0;
    dfdtheta(7,2) = 0;
    
end