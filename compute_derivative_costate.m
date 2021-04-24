function lambda_dot = compute_derivative_costate(lambda,Z, theta, noise_on, fixed_prob)
%COMPUTE_DERIVATIVE_COSTATE Summary of this function goes here
%   Detailed explanation goes here
dfdx = compute_dfdx_system(Z,theta, noise_on, fixed_prob);
lambda_dot = -(lambda'*dfdx)';
end

