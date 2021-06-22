function zdot = compute_derivative_system(Z, input, noise_on, fixed_prob, pars)
    
    R = Z(1);
    S = Z(2);
    E = Z(3);
    L = Z(4);
    I = Z(5);
    V = Z(6);
    A = Z(7);

    %% compute monod and probability
    Psi = monod(R,pars.mu_max,pars.R_in);
%     for i=1:sample_size
%         P_A = probability(noisy_A(A),input,A0);
%         sum = sum+P_A;
%     end    
%     fraction_of_E = sum*E/sample_size;     
    P_A_cumulative = probability(A,input,pars.A0, noise_on, fixed_prob);
    
    %% diff eq
    R_dot = pars.J - pars.e*Psi*(L+(1-pars.alpha_S)*S)-pars.d_R*R;
    S_dot = Psi*(1-pars.alpha_S)*S-pars.phi*S*V-pars.d_S*S;
    E_dot = pars.phi*S*V - pars.lambda*E - pars.d_E*E;
    L_dot = Psi*L + P_A_cumulative*E*pars.lambda - pars.d_L*L;
    I_dot = (1-P_A_cumulative)*E*pars.lambda - pars.eta*I - pars.d_I*I;
    V_dot = pars.beta*pars.eta*I - pars.phi*S*V - pars.d_V*V - pars.phi*(E+L+I)*V;
    A_dot = pars.k_A_I*I + pars.k_A_L*L - pars.d_A*A;
    
    zdot = [R_dot,S_dot,E_dot,L_dot,I_dot,V_dot,A_dot]';
end
