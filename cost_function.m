function J = cost_function(Z0, u, t, dt, noise_on, fixed_prob)
    max_iter = 1;
    J = 0;
    for i = 1:max_iter
        Z = forward_euler(Z0, u, t, dt, noise_on, fixed_prob);
        J = J + (Z(4,end)+Z(3,end));
    end
    J = J/max_iter;
end