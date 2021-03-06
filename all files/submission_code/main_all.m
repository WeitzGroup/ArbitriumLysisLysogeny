close all

option = 3; % set 1 for iterating over r0 and J, 2 for iterating over production rates, 3 for low production rate
if option == 1 | option == 3
    prod_rate_vector = 5;  % fixed production rate
    tf_vector = [12,18,24,30,36,42,48,72];
elseif option == 2    
    prod_rate_vector = 5:0.25:8; % to sweep for different generation rates
    tf_vector = 24;
end

if option == 2 | option ==3
    J_val_vec = 0;
    r0_vec = 40;
else
    J_val_vec = [0, 1, 2, 3];
    r0_vec = [40, 60, 80, 100];
end

folder=pwd;
filepath = strcat(folder,'/data__files');
addpath(genpath(strcat(folder,'/helper_functions')));

optimal_params_save = zeros(length(prod_rate_vector),2,length(tf_vector));
optimal_params_no_production = zeros(2, length(tf_vector));

noise_on = 0;
fixed_prob = 0;

%% params
% known values
for J_val = J_val_vec    
for r0 = r0_vec  
pars.mu_max = 1.2;
pars.R_in = 4;
pars.J = J_val;  
pars.e = 5*10^(-7); 
pars.d_R = 0.0;  
pars.d_S = 0.075;
pars.d_E = 0.075;
pars.d_L = 0.075;
pars.d_I = 0.075;
pars.lambda = 2;
pars.eta = 1;
pars.phi = 3.4*10^(-10);
pars.alpha_S = -0.1;

pars.d_V = 0.075;              
pars.d_A = 0.1;
pars.k_A_L = 5*10^7; 
pars.k_A_I = 5*10^7;
pars.A0 = 10^12; 
pars.beta = 50;

if option == 3
    pars.k_A_L = 10^7; 
    pars.k_A_I = 10^7; 
end    

prod_rate_iter = 0;    
for prod_rate = prod_rate_vector    %(loop over different production rates)
prod_rate_iter = prod_rate_iter + 1;    
%disp(prod_rate)

if option ==2
    pars.k_A_L = 10^(prod_rate); %(only for sweeping over different production rates)
    pars.k_A_I = 10^(prod_rate);
end

for tf_index = 1:length(tf_vector) %(loop over different time horizons)
%disp(tf_index)
%% time vector
dt = 20/3600;    % dt = 20s, tf = 12 hours
t = 0:dt:tf_vector(tf_index);

%% initial concentrations
R0 = r0;
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
J_brute = zeros(length(sigmoid_shaping_values), length(switch_point_values));
for sigmoid_shaping = sigmoid_shaping_values
    index_1 = index_1+1;
    for switch_point = switch_point_values
        index_2 = index_2 +1;
        J_brute(index_1,index_2) = cost_function(Z0, [sigmoid_shaping;switch_point], t, dt, noise_on, fixed_prob, pars);
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
J = 10^10*ones(1, 30);
J(loop_iter) = cost_function(Z0, optimal_input_old, t, dt, noise_on, fixed_prob, pars);

%% termination condition
change_in_J = J(loop_iter);
cut_off_for_J_change = 10^(-5);
while change_in_J > cut_off_for_J_change && loop_iter < 30 % termination condition
    loop_iter = loop_iter+1;
    inputs_iterated(:,loop_iter) = optimal_input_old;

    %% step 1
    Z = forward_euler(Z0, optimal_input_old, t, dt, noise_on, fixed_prob, pars);            % state trajectory (fwd)
    
    %% step 2
    lambda_tf = [0,0,1,1,0,0,0]';
    lambda = backward_euler_costate(lambda_tf, Z, optimal_input_old, t, dt, noise_on, fixed_prob, pars);
    
    %% step 3
    dphi_dx = [0,0,0,1,1,0,0];
    grad_j = grad_J(lambda, Z, optimal_input_old, t, dt, noise_on, fixed_prob, pars);
%     if norm(grad_j)>10
%         grad_j = 10*grad_j / norm(grad_j);             % normalize?
%     end
    %% step 4
    new_input = minimum_step(Z0, grad_j, optimal_input_old,t,dt,alpha,beta, noise_on, fixed_prob, pars);
    
    J(loop_iter) = cost_function(Z0, new_input, t, dt, noise_on, fixed_prob, pars);
    change_in_J = J(loop_iter) - J(loop_iter-1);    % update termination condition
    if change_in_J>0
        optimal_input_old = new_input;                  % assign to optimal_input_old for next iteration
        J_max = J(loop_iter);
    end
end    
optimal_params_save(prod_rate_iter,:,tf_index) = optimal_input_old;
optimal_params_no_production(:,tf_index) = optimal_params_save(1,:,tf_index);

end

optimal_params = optimal_params_no_production;


if option == 1
    name = ['optimal_params_r0_' num2str(R0) '_j_' num2str(pars.J) '.mat'];
elseif option ==2
    name = 'optimal_params_r0_40_j_0_sweep.mat';   % for different production rates
elseif option ==3
    name = 'optimal_params_r0_40_j_0_low_production.mat'; % for low production rate
end

folder=pwd;
filepath = strcat(folder,'/data__files');
matname = fullfile(filepath, name);
if option == 1 | option ==3
    save(matname, 'optimal_params')
elseif option ==2
    save(matname, 'optimal_params_save')
end

end


end 
end