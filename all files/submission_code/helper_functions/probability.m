function P_A = probability(A_true, optimal_input, switch_pt, noise_on, fixed_prob)
if fixed_prob
    if fixed_prob == 0.01
        P_A = 0;
    else    
        P_A = fixed_prob;
    end    
else    
    P_A = 0.0;
    max_val = 1;
    A0 = switch_pt;
    k = optimal_input(1);
    x0 = optimal_input(2);
%    P_A = optimal_input;
    lambda = A_true/10^12;
    if noise_on
        sigma = 1/(1+exp(-k*(A_true/(A0*x0)-1)));
        P_A = max_val/(1+exp(-k*(A_true/(A0*x0)-1))) + (A_true/2)*sigma*(1-sigma)*(1-2*sigma)*max_val*(k/(A0*x0))^2; 
    else
        P_A = max_val/(1+exp(-k*(A_true/(A0*x0)-1))); 
    end 
end    
end