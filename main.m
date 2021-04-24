close all
clear all

%% time vector
dt = 20/3600; tf =48;   % dt = 20s, tf = 12 hours
t = 0:dt:tf;

%% initial concentrations
R0 = 10;
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

sigmoid_shaping_values = [-10 0.0001 0.001 0.01 0.1 0.5 1 2 10];
switch_point_values = [10^(-8) 0.0001 0.001 0.01 0.1 0.5 1 5 10 30 10^5];
%optimal_input_brute = [1;1];
for sigmoid_shaping = sigmoid_shaping_values
    index_1 = index_1+1
    for switch_point = switch_point_values
        index_2 = index_2 +1
        J_brute(index_1,index_2) = cost_function(Z0, [sigmoid_shaping;switch_point], t, dt);
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
J(loop_iter) = cost_function(Z0, optimal_input_old, t, dt);

%% termination condition
change_in_J = J(loop_iter);
cut_off_for_J_change = 10^(-5);
while change_in_J > cut_off_for_J_change && loop_iter < 30
    loop_iter = loop_iter+1
    inputs_iterated(:,loop_iter) = optimal_input_old;

    %% step 1
    disp('step 1')
    Z = forward_euler(Z0, optimal_input_old, t, dt);            % state trajectory (fwd)
    
    %% step 2
    disp('step 2')
    lambda_tf = [0,0,1,1,0,0,0]';
    lambda = backward_euler_costate(lambda_tf, Z, optimal_input_old, t, dt);
    
    %% step 3
    disp('step 3')
    dphi_dx = [0,0,0,1,1,0,0];
    grad_j = grad_J(lambda, Z, optimal_input_old, t, dt);
    if norm(grad_j)>1
        grad_j = grad_j / norm(grad_j);             % normalize?
    end
    %% step 4
    disp('step 4')    
    new_input = minimum_step(Z0, grad_j, optimal_input_old,t,dt,alpha,beta);
    
    J(loop_iter) = cost_function(Z0, new_input, t, dt);
    change_in_J = J(loop_iter) - J(loop_iter-1);    % update termination condition
    if change_in_J>0
        optimal_input_old = new_input;                  % assign to optimal_input_old for next iteration
        J_max = J(loop_iter);
    end
end    



[param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14,param15,param16,param17,param18] = return_parameters(optimal_input_old);
probability_switch_pt_initial_guess = param18;
k = optimal_input_old(1);
x0 = optimal_input_old(2);

Z = forward_euler(Z0, optimal_input_old, t, dt);



for i=1:length(t)
    prob(i) =  probability(Z(7,i), optimal_input_old, probability_switch_pt_initial_guess);
end

line_thickness = 2;

figure(999)
plot(t, prob, 'Linewidth', line_thickness)
xlabel('$Time~[h]$','Interpreter','latex')
ylabel('$Probability$','Interpreter','latex')
fig.PaperUnits = 'inches';
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
fig.PaperUnits = 'inches';
pbaspect([2.5 1 1])

figure(2)
for i = 0:50
    span_A(i+1) = 10^(6 + i/50*8);
end
prob_A = zeros(1,length(span_A));
for i = 1:length(span_A)
    prob_A(i) = probability(span_A(i), optimal_input_old, probability_switch_pt_initial_guess);
end    
plot (log10(span_A),prob_A, 'Linewidth', line_thickness)
%title ('Probability as a fn of A')
xlabel('$A~[molecules/ml]$','Interpreter','latex')
ylabel('$Probability$','Interpreter','latex')
ylim ([0,1]);
xticks([6 7 8 9 10 11 12 13 14])
xticklabels({'10^{6}','10^{7}','10^{8}','10^{9}','10^{10}','10^{11}','10^{12}','10^{13}','10^{14}'})
%xticks([10 10.5 11 11.5 12 12.5 13 13.5 14])
%xticklabels({'10^{10}','10^{10.5}','10^{11}','10^{11.5}','10^{12}','10^{12.5}','10^{13}','10^{13.5}','10^{14}'})
fig.PaperUnits = 'inches';
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
fig.PaperUnits = 'inches';
pbaspect([2.5 1 1])




% Z = forward_euler(Z0, optimal_input_old, t, dt);
figure (1)
%plot(t,log10(Z(1,:)))
hold on
plot(t,log10(Z(2,:)), 'Linewidth', line_thickness)
%plot(t,log10(Z(3,:)))
plot (t,log10(Z(4,:)), 'Linewidth', line_thickness)
plot(t,log10(Z(5,:)), 'Linewidth', line_thickness)
plot(t,log10(Z(6,:)), 'Linewidth', line_thickness)
hold off
%legend('R','S','E','L','I','V','A')
legend('S','L','I', 'V')
ylim([1,9]);
yticks([1 2 3 4 5 6 7 8 9])
yticklabels({'10^1','10^2','10^3','10^4','10^5','10^6','10^7','10^8','10^9'})
xlabel('$Time~[h]$','Interpreter','latex')
ylabel('$Population~[ml^{-1}]$','Interpreter','latex')
fig.PaperUnits = 'inches';
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
fig.PaperUnits = 'inches';
pbaspect([2.5 1 1])





