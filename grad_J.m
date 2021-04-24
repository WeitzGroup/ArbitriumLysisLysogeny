function grad_j= grad_J(lambda, Z, optimal_input_old, t, dt, noise_on, fixed_prob)
    grad_j = zeros(1,2);
    
    for i = 1:length(t)
        dfdtheta = compute_dfdtheta(Z(:,i),optimal_input_old, noise_on, fixed_prob);
        grad_j = grad_j + lambda(:,i)'*dfdtheta*dt;
    end    
end

