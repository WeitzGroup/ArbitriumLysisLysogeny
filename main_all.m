close all

%tf_vector = [2,6,12,24,36,48,72,96,120,240];
%tf_vector = [4,8,10,12,14,16];
tf_vector = [12,18,24,30,36,42,48,72];
%load('optimal_params_sensing_noise_latest.mat')
optimal_params_save = zeros(2,length(tf_vector));
noise_on = 0;
fixed_prob = 0;
for tf_index = 1:length(tf_vector)
disp(tf_index)
%% time vector
dt = 20/3600;    % dt = 20s, tf = 12 hours
t = 0:dt:tf_vector(tf_index);

%% initial concentrations
R0 = 40;%25 
S0 = 10^7;
E0 = 0;
L0 = 0;
I0 = 0;
V0 = 10^5;
A0 = 0;

Z0 = [R0,S0,E0,L0,I0,V0,A0];

%% vector to store states 
Z = zeros(7,length(t));

%% costate vector and boundary condition
P = zeros(7, length(t));

%% optimal input
index_1 = 0; index_2 = 0;
J_max = 0;

sigmoid_shaping_values = [-10 0.0001 1 10 ];
switch_point_values = [10^(-8) 0.0001 0.1 0.3 1 2 5 10 20 30 50 100 200 1000];
%optimal_input_brute = [1;1];
for sigmoid_shaping = sigmoid_shaping_values
    index_1 = index_1+1;
    for switch_point = switch_point_values
        index_2 = index_2 +1;
        J_brute(index_1,index_2) = cost_function(Z0, [sigmoid_shaping;switch_point], t, dt, noise_on, fixed_prob);
        if J_brute(index_1,index_2)>J_max
            J_max = J_brute(index_1,index_2);
            optimal_input_brute = [sigmoid_shaping;switch_point];
        end
    end
    index_2 = 0;
end

optimal_input_old = optimal_input_brute;      % k = 1 and x0 = 150 nM and k_A = 5*10^9 molecules/(cell*ml)
inputs_iterated = zeros(2,100);

%% set constants for Armijo Algo
alpha = 10^(-3);
beta = 0.3;

%% N iterations for converging to optimal input (change to normed dist b/w u and u_hat less than delta later)
N = 5;
cut_off_for_input_change = 0.01;
loop_iter = 1;
J(loop_iter) = cost_function(Z0, optimal_input_old, t, dt, noise_on, fixed_prob);

%% termination condition
change_in_J = J(loop_iter);
cut_off_for_J_change = 10^(-5);
while change_in_J > cut_off_for_J_change && loop_iter < 30
    loop_iter = loop_iter+1;
    inputs_iterated(:,loop_iter) = optimal_input_old;

    %% step 1
    Z = forward_euler(Z0, optimal_input_old, t, dt, noise_on, fixed_prob);            % state trajectory (fwd)
    
    %% step 2
    lambda_tf = [0,0,1,1,0,0,0]';
    lambda = backward_euler_costate(lambda_tf, Z, optimal_input_old, t, dt, noise_on, fixed_prob);
    
    %% step 3
    dphi_dx = [0,0,0,1,1,0,0];
    grad_j = grad_J(lambda, Z, optimal_input_old, t, dt, noise_on, fixed_prob);
%     if norm(grad_j)>10
%         grad_j = 10*grad_j / norm(grad_j);             % normalize?
%     end
    %% step 4
    new_input = minimum_step(Z0, grad_j, optimal_input_old,t,dt,alpha,beta, noise_on, fixed_prob);
    
    J(loop_iter) = cost_function(Z0, new_input, t, dt, noise_on, fixed_prob);
    change_in_J = J(loop_iter) - J(loop_iter-1);    % update termination condition
    if change_in_J>0
        optimal_input_old = new_input;                  % assign to optimal_input_old for next iteration
        J_max = J(loop_iter);
    end
end    
optimal_params_save(:,tf_index) = optimal_input_old;
%save('optimal_params_no_noise_latest.mat','optimal_params_save')
save('optimal_params_r0_40_j_3.mat','optimal_params_save')
end



