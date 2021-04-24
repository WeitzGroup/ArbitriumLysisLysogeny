function dfdx = compute_dfdx_system(x_hat, u_hat,noise_on, fixed_prob)
    
    %% Assign letters from x_hat
    R = x_hat(1);
    S = x_hat(2);
    E = x_hat(3);
    L = x_hat(4);
    I = x_hat(5);
    V = x_hat(6);
    A = x_hat(7);

    %% parameters assigned
    [mu_max,R_in,J,e,d_R,d_S,d_E,d_L,d_I,lambda,eta,beta,phi,alpha_S,d_V,d_A,k_A_mean,A0] = return_parameters(u_hat);
    k_A = k_A_mean;
    k_L = k_A;
    %% compute monod and probability
    Psi = monod(R,mu_max,R_in);
    P_A = probability(A,u_hat,A0,noise_on, fixed_prob);
    dfdx = zeros(7,7);
    
    %% compute needed derivatives
    dPsidR = mu_max/(R+R_in) - mu_max*R/((R+R_in)^2);
    
    max_val = 1; k = u_hat(1); x0 = u_hat(2);
    logistic_fn = 1/(1+exp(-k*(A/(A0*x0)-1)));
    dP_AdA = max_val*k/(A0*x0)*(logistic_fn*(1-logistic_fn));
    
    %% Assign value to dfdx term by term
    
    dfdx(1,1) = - e*dPsidR*(L+(1-alpha_S)*S) - d_R;
    dfdx(1,2) = - e*Psi*(1-alpha_S);
    dfdx(1,3) = 0;
    dfdx(1,4) = - e*Psi;
    dfdx(1,5) = 0;
    dfdx(1,6) = 0;
    dfdx(1,7) = 0;
    
    dfdx(2,1) = dPsidR*(1-alpha_S)*S;
    dfdx(2,2) = Psi*(1-alpha_S) - phi*V - d_S;
    dfdx(2,3) = 0;
    dfdx(2,4) = 0;
    dfdx(2,5) = 0;
    dfdx(2,6) = - phi*S;
    dfdx(2,7) = 0;
    
    dfdx(3,1) = 0;
    dfdx(3,2) = phi*V;
    dfdx(3,3) = - lambda - d_E;
    dfdx(3,4) = 0;
    dfdx(3,5) = 0;
    dfdx(3,6) = phi*S;
    dfdx(3,7) = 0;
    
    dfdx(4,1) = dPsidR*L;
    dfdx(4,2) = 0;
    dfdx(4,3) = P_A*lambda;
    dfdx(4,4) = Psi - d_L;
    dfdx(4,5) = 0;
    dfdx(4,6) = 0;
    dfdx(4,7) = dP_AdA*lambda*E;
    
    dfdx(5,1) = 0;
    dfdx(5,2) = 0;
    dfdx(5,3) = (1-P_A)*lambda ;
    dfdx(5,4) = 0;
    dfdx(5,5) = - eta - d_I;
    dfdx(5,6) = 0;
    dfdx(5,7) = -dP_AdA*lambda*E;
    
    dfdx(6,1) = 0;
    dfdx(6,2) = -phi*V;
    dfdx(6,3) = -phi*V;
    dfdx(6,4) = -phi*V;
    dfdx(6,5) = beta*eta - phi*V;
    dfdx(6,6) = - phi*S - d_V;
    dfdx(6,7) = 0;
    
    dfdx(7,1) = 0;
    dfdx(7,2) = 0;
    dfdx(7,3) = 0;
    dfdx(7,4) = k_L;
    dfdx(7,5) = k_A;
    dfdx(7,6) = 0;
    dfdx(7,7) = - d_A;
end