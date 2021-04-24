function zdot = compute_derivative_system(Z, input, noise_on, fixed_prob)
    
    R = Z(1);
    S = Z(2);
    E = Z(3);
    L = Z(4);
    I = Z(5);
    V = Z(6);
    A = Z(7);
    %% parameters assigned
    [mu_max,R_in,J,e,d_R,d_S,d_E,d_L,d_I,lambda,eta,beta,phi,alpha_S,d_V,d_A,k_A_mean,A0] = return_parameters(input);   
    
    k_A = k_A_mean;
    k_L = k_A;
    %% compute monod and probability
    Psi = monod(R,mu_max,R_in);
%     for i=1:sample_size
%         P_A = probability(noisy_A(A),input,A0);
%         sum = sum+P_A;
%     end    
%     fraction_of_E = sum*E/sample_size;     
    P_A_cumulative = probability(A,input,A0, noise_on, fixed_prob);
    
    %% diff eq
    R_dot = J - e*Psi*(L+(1-alpha_S)*S)-d_R*R;
    S_dot = Psi*(1-alpha_S)*S-phi*S*V-d_S*S;
    E_dot = phi*S*V - lambda*E - d_E*E;
    L_dot = Psi*L + P_A_cumulative*E*lambda - d_L*L;
    I_dot = (1-P_A_cumulative)*E*lambda - eta*I - d_I*I;
    V_dot = beta*eta*I - phi*S*V - d_V*V - phi*(E+L+I)*V;
    A_dot = k_A*I + k_L*L - d_A*A;
    
    zdot = [R_dot,S_dot,E_dot,L_dot,I_dot,V_dot,A_dot]';
end
