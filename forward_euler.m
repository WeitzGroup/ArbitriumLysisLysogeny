%% function to compute state trajectory given an input
function Z = forward_euler(Z0, optimal_input, t, dt, noise_on, fixed_prob)
    Z = zeros(7,length(t));
    Z(:,1) = Z0;

    %% run forward euler
    for i=2:length(t)
        Zdot = compute_derivative_system(Z(:,i-1),optimal_input, noise_on, fixed_prob);
        Z(:,i)= Z(:,i-1) + Zdot*dt;           % change per hour to per second
    end    
end



