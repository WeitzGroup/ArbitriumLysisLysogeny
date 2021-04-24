function xidot = compute_derivative_xi(xi,Z, theta, noise_on, fixed_prob)
    dfdx = compute_dfdx_system(Z,theta,noise_on, fixed_prob);
    dfdtheta = compute_dfdtheta(Z,theta, noise_on, fixed_prob);
    xidot = dfdx * xi + dfdtheta;
end
