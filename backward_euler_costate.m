function lambda = backward_euler_costate(lambda_tf,Z, optimal_input, t, dt, noise_on, fixed_prob)
    lambda = zeros(7,length(t));
    lambda(:,end) = lambda_tf;
    
    %% run forward euler
    for i=length(t)-1:-1:1
        lambda_dot = compute_derivative_costate(lambda(:,i+1),Z(:,i+1),optimal_input, noise_on, fixed_prob);
        lambda(:,i)= lambda(:,i+1) - lambda_dot*dt;           
    end    

end

