function u_new = minimum_step(Z0, grad_J, u, t, dt, alpha, beta, noise_on, fixed_prob)

    J_old = cost_function(Z0, u, t, dt, noise_on,  fixed_prob);
    diff_max = 0.0;
    idx_max = -1;
%    while J_new - J_old <= alpha*beta^k*norm(grad_J)
    for k = 0:40
        J_new = cost_function(Z0, u + beta^k*grad_J', t, dt, noise_on, fixed_prob);
        
        diff = J_new - J_old;
        
        if (diff>diff_max)
            diff_max = diff;
            idx_max = k;
        end    
    end 
    
%     if k>0
%         k = k-1;        % new value of k from while loop doesn't satisfy condition, so take prev value
%     end
    if idx_max>=0
        u_new = u + beta^idx_max*grad_J';
    else
        u_new = u;
    end
end
