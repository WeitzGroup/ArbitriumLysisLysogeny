function dfdx = compute_dfdx_system(x_hat, u_hat,noise_on, fixed_prob, pars)
    
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
    P_A = probability(A,u_hat,pars.A0,noise_on, fixed_prob);
    dfdx = zeros(7,7);
    
    %% compute needed derivatives
    dPsidR = pars.mu_max/(R+pars.R_in) - pars.mu_max*R/((R+pars.R_in)^2);
    
    max_val = 1; k = u_hat(1); x0 = u_hat(2);
    logistic_fn = 1/(1+exp(-k*(A/(pars.A0*x0)-1)));
    dP_AdA = max_val*k/(pars.A0*x0)*(logistic_fn*(1-logistic_fn));
    
    %% Assign value to dfdx term by term
    
    dfdx(1,1) = - pars.e*dPsidR*(L+(1-pars.alpha_S)*S) - pars.d_R;
    dfdx(1,2) = - pars.e*Psi*(1-pars.alpha_S);
    dfdx(1,3) = 0;
    dfdx(1,4) = - pars.e*Psi;
    dfdx(1,5) = 0;
    dfdx(1,6) = 0;
    dfdx(1,7) = 0;
    
    dfdx(2,1) = dPsidR*(1-pars.alpha_S)*S;
    dfdx(2,2) = Psi*(1-pars.alpha_S) - pars.phi*V - pars.d_S;
    dfdx(2,3) = 0;
    dfdx(2,4) = 0;
    dfdx(2,5) = 0;
    dfdx(2,6) = - pars.phi*S;
    dfdx(2,7) = 0;
    
    dfdx(3,1) = 0;
    dfdx(3,2) = pars.phi*V;
    dfdx(3,3) = - pars.lambda - pars.d_E;
    dfdx(3,4) = 0;
    dfdx(3,5) = 0;
    dfdx(3,6) = pars.phi*S;
    dfdx(3,7) = 0;
    
    dfdx(4,1) = dPsidR*L;
    dfdx(4,2) = 0;
    dfdx(4,3) = P_A*pars.lambda;
    dfdx(4,4) = Psi - pars.d_L;
    dfdx(4,5) = 0;
    dfdx(4,6) = 0;
    dfdx(4,7) = dP_AdA*pars.lambda*E;
    
    dfdx(5,1) = 0;
    dfdx(5,2) = 0;
    dfdx(5,3) = (1-P_A)*pars.lambda ;
    dfdx(5,4) = 0;
    dfdx(5,5) = - pars.eta - pars.d_I;
    dfdx(5,6) = 0;
    dfdx(5,7) = -dP_AdA*pars.lambda*E;
    
    dfdx(6,1) = 0;
    dfdx(6,2) = -pars.phi*V;
    dfdx(6,3) = -pars.phi*V;
    dfdx(6,4) = -pars.phi*V;
    dfdx(6,5) = pars.beta*pars.eta - pars.phi*V;
    dfdx(6,6) = - pars.phi*S - pars.d_V;
    dfdx(6,7) = 0;
    
    dfdx(7,1) = 0;
    dfdx(7,2) = 0;
    dfdx(7,3) = 0;
    dfdx(7,4) = pars.k_A_L;
    dfdx(7,5) = pars.k_A_I;
    dfdx(7,6) = 0;
    dfdx(7,7) = - pars.d_A;
end