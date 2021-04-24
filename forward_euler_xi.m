function xi = forward_euler_xi(xi_0, Z, optimal_input, t, dt)
    xi = zeros(7,2,length(t));
    xi(:,:,1) = xi_0;

    %% run forward euler
    for i=2:length(t)
        xi_dot = compute_derivative_xi(xi(:,:,i-1),Z(:,i-1),optimal_input);
        xi(:,:,i)= xi(:,:,i-1) + xi_dot*dt;           
    end    
end

