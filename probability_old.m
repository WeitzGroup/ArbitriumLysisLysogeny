function P_A = probability_old(A_true, optimal_input, switch_pt, noise_on, fixed_prob)
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
        if lambda < 5
            for i = 0:30
                prob = lambda^i*exp(-lambda)/factorial(i);
                P_A = prob*max_val/(1+exp(-k*(i*10^12/(A0*x0)-1))) + P_A;
            end
        elseif lambda >= 5 
            for i = 0:6*floor(lambda)
                if P_A >= 0.99999999 
                    P_A = 1;
                    break
                end
                prob = exp(-lambda/2);
                for j = 1:i
                    prob = prob*lambda/j;   
                end    
                prob = prob*exp(-lambda/2);
                P_A = prob*max_val/(1+exp(-k*(i*10^12/(A0*x0)-1))) + P_A;            
            end
        else
            
        end
    else
        P_A = max_val/(1+exp(-k*(A_true/(A0*x0)-1))); 
    end 
end    
end