close all

noise_on = 1;

% get current path
folder=pwd;
% Addpath to subfolders
addpath(genpath(strcat(folder,'/data__files')));


load('optimal_params_r0_40_j_0.mat')
%optimal_params = optimal_params_no_production;

%% params
%% params
% known values
pars.mu_max = 1.2;%1.2;
pars.R_in = 4;%4;
pars.J = 0;  % 0.2 prev
pars.e = 5*10^(-7); %5*10^(-7);
pars.d_R = 0.0;  % 0.2 prev
pars.d_S = 0.075;
pars.d_E = 0.075;
pars.d_L = 0.075;
pars.d_I = 0.075;
pars.lambda = 2;
pars.eta = 1;
pars.phi = 3.4*10^(-10);
pars.alpha_S = -0.1;

% doubtful
pars.d_V = 0.075;              % table says 0.4?
% unknown
pars.d_A = 0.1;
pars.k_A_L = 5*10^7; % 6.022*10^6
pars.k_A_I = 5*10^7;
pars.A0 = 10^12; 

% beta = tradeoff_value(input(3));
pars.beta = 50;

probability_switch_pt_initial_guess = pars.A0;
%tf_vector = [2,6,12,24,36,48,72,96,120,240];
%tf_vector = [4,8,10,12,14,16];
tf_vector = [12,18,24,30,36,42,48,72];

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%color_cell = {'k','k','k','k'};
ls = {'-','--',':','-.'};
ms = {'o','s','d','h'};
line_thickness = 3;

%% initial concentrations
R0 = 40;
S0 = 10^7;
E0 = 0;
L0 = 0;
I0 = 0;
V0 = 10^5;
%V0 = 0;
A0 = 0;

Z0 = [R0,S0,E0,L0,I0,V0,A0];
dt = 20/3600; 

%% figure 2(a), 2(b) popuation dynamics for fixed strats for 48 hours
mixed_prob = 0.5;   
t = 0:dt:48;
index_for_48_hours = 7;
Z = forward_euler(Z0, optimal_params(:,index_for_48_hours), t, dt, noise_on, 0, pars);
Z_48 = Z;
%save('Z_48.mat','Z')    
Z_fixed_lysis = forward_euler(Z0, optimal_params(:,index_for_48_hours), t, dt, noise_on, 0.01, pars);
Z_fixed_mixed = forward_euler(Z0, optimal_params(:,index_for_48_hours), t, dt, noise_on, mixed_prob, pars);
Z_fixed_lysogeny = forward_euler(Z0, optimal_params(:,index_for_48_hours), t, dt, noise_on, 1, pars);

% lysis only
figure;
subplot(1,3,1)
plot(t, (Z_fixed_lysis(2,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_lysis(3,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_lysis(4,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_lysis(6,:)),  'Linewidth', line_thickness)
hold on
%plot(t, log10(Z_fixed_lysis(7,:)),  'Linewidth', line_thickness)    
%xlabel('$Time~[hours]$','Interpreter','latex')
ylabel('$\mathrm{Population~[ml^{-1}}]$','Interpreter','latex')
xticks([0,24,48]);
ylim([1e3, 1e9]);
yticks([1e3  1e5  1e7 1e9])
xlim([0,48])
title('P = 0','Interpreter','latex')
%legend('S','E','L','V')
set(gca,'FontSize',20);
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')


fig.PaperUnits = 'inches';
pbaspect([2.5 2.5 1])

% lysogeny only
subplot(1,3,3)
plot(t, (Z_fixed_lysogeny(2,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_lysogeny(3,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_lysogeny(4,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_lysogeny(6,:)),  'Linewidth', line_thickness)
hold on
%plot(t, log10(Z_fixed_lysogeny(7,:)),  'Linewidth', line_thickness)
%xlabel('$Time~[hours]$','Interpreter','latex')
xticks([0,24,48]);
xlim([0,48])
ylim([1e3, 1e9]);
yticks([1e3  1e5  1e7 1e9])
yticklabels({'10^3','10^5','10^7','10^9'})    
%ylabel('$Population~[ml^{-1}]$','Interpreter','latex')
title('P = 1','Interpreter','latex')
legend('S','E','L','V','Interpreter','latex')
legend boxoff
set(gca,'FontSize',20);
fig.PaperUnits = 'inches';
pbaspect([2.5 2.5 1])
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')

% mixed strat
subplot(1,3,2)
plot(t, (Z_fixed_mixed(2,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_mixed(3,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_mixed(4,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_mixed(6,:)),  'Linewidth', line_thickness)
hold on
%plot(t, log10(Z_fixed_half(7,:)),  'Linewidth', line_thickness)
xlabel('Time [hours]','Interpreter','latex')
xticks([0,24,48]);
xlim([0,48])
ylim([1e3, 1e9]);
yticks([1e3  1e5  1e7 1e9])
yticklabels({'10^3','10^5','10^7','10^9'})    
%ylabel('$Population~[ml^{-1}]$','Interpreter','latex')
title('P = 0.5','Interpreter','latex')
%legend('S','E','L','V')    
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gca,'FontSize',20);
fig.PaperUnits = 'inches';
pbaspect([2.5 2.5 1])
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')

% compare
figure;
%plot(t, log10(Z(3,:)+Z(4,:)),  'Linewidth', line_thickness)
%hold on
plot(t, (Z_fixed_lysis(3,:)+Z_fixed_lysis(4,:)), 'k--', 'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_mixed(3,:)+Z_fixed_mixed(4,:)),'k:',  'Linewidth', line_thickness)
hold on
plot(t, (Z_fixed_lysogeny(3,:)+Z_fixed_lysogeny(4,:)),'k-.',  'Linewidth', line_thickness)

fig.PaperUnits = 'inches';
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gca,'FontSize',20);
%set(gca,'TickLabelInterpreter', 'latex');
fig.PaperUnits = 'inches';
pbaspect([2.5 1 1])
xlabel('Time [hours]','Interpreter','latex')
xticks([0,6,12,18,24,30,36,42,48]);
xlim([12,48])
ylim([1e2, 1e6]);
yticks([1e2 1e3  1e4  1e5  1e6 ])

ylabel('$\mathrm{Birth~State~Population~[ml^{-1}]}$','Interpreter','latex')
%title('$Total~birth~states~comparison$','Interpreter','latex')
legend('P = 0','P = 0.5','P = 1','Interpreter','latex')    
pbaspect([2.5 1.5 1])
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')
legend boxoff


%% figure 3(a) in paper; prob vs A for different tf

span_A = zeros(1, 1001);
for i = 0:1000
    span_A(i+1) = 10^(8 + i/1000*8);
end

conversion_factor = 10^3*10^9/(6.022*10^23);    % ml, nano, mole
span_A_nanomolar = span_A*conversion_factor;

figure;
prob_A = zeros(length(optimal_params),length(span_A));
prob_final_times = [12,18,24,36,48];
%color_cell = {'r','y','g','','c','',[250/255,107/255,242/255],''};
color_cell = {[0,0,0],[0.2,0.2,0.2],[0.4,0.4,0.4],'',[0.6,0.6,0.6],'',[0.8,0.8,0.8],''};
for tf_idx = 1:length(optimal_params)
    if (ismember(tf_vector(tf_idx), prob_final_times))   
    for i = 1:length(span_A)
        prob_A(tf_idx, i) = probability(span_A(i), optimal_params(:,tf_idx), probability_switch_pt_initial_guess, noise_on, 0);
    end   
    
    plot (span_A_nanomolar, prob_A(tf_idx,:), 'color', color_cell{tf_idx},'Linewidth', line_thickness)
    hold on
    end
end    

xlabel('Arbitrium~concentration, $A$ ~$\mathrm{[nM]}$','Interpreter','latex')
ylabel('Probability~of~lysogeny, $P_{opt}(A)$','Interpreter','latex')
ylim ([0,1]);
xlim([1e0,1e3]);
xticks([ 1 10 100 1000])
xticklabels({1,10,100,1000})


set(gca,'TickLabelInterpreter','latex')
fig.PaperUnits = 'inches';
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gca,'FontSize', 30);
set(gca,'xscale','log')

% set(gca,'TickLabelInterpreter', 'latex');
fig.PaperUnits = 'inches';
pbaspect([2.5 1.5 1])
legend('12 hours','18 hours','24 hours','36 hours','48 hours','Interpreter','latex')
legend boxoff

%% figure 3(b) in paper; prob vs A for different tf (low production)
load('optimal_params_r0_40_j_0_low_production.mat')

figure;
prob_A = zeros(length(optimal_params),length(span_A));
prob_final_times = [12,18,24,36,48];
%color_cell = {'r','y','g','','c','',[250/255,107/255,242/255],''};
color_cell = {[0,0,0],[0.2,0.2,0.2],[0.4,0.4,0.4],'',[0.6,0.6,0.6],'',[0.8,0.8,0.8],''};
for tf_idx = 1:length(optimal_params)
    if (ismember(tf_vector(tf_idx), prob_final_times))   
    for i = 1:length(span_A)
        prob_A(tf_idx, i) = probability(span_A(i), optimal_params(:,tf_idx), probability_switch_pt_initial_guess, noise_on, 0);
    end   
    
    plot (span_A_nanomolar, prob_A(tf_idx,:), 'color', color_cell{tf_idx},'Linewidth', line_thickness)
    hold on
    end
end    

xlabel('Arbitrium~concentration, $A$ ~$\mathrm{[nM]}$','Interpreter','latex')
ylabel('Probability~of~lysogeny, $P_{opt}(A)$','Interpreter','latex')
ylim ([0,1]);

xlim([1e-1,1e2]);         % use with low production rate
xticks([1e-1 1 10 100])
xticklabels({0.1,1,10,100})

set(gca,'TickLabelInterpreter','latex')
fig.PaperUnits = 'inches';
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gca,'FontSize', 30);
set(gca,'xscale','log')

% set(gca,'TickLabelInterpreter', 'latex');
fig.PaperUnits = 'inches';
pbaspect([2.5 1.5 1])
legend('12 hours','18 hours','24 hours','36 hours','48 hours','Interpreter','latex')
legend boxoff

load('optimal_params_r0_40_j_0.mat')

%% Fig. 4 population dynamics for different tf (12, 18, 24, 48)

Z_temp = zeros(length(tf_vector),7,length(0:dt:48));
for i = 1:length(tf_vector)
    t = 0:dt:tf_vector(i);
    Z = forward_euler(Z0, optimal_params(:,i), t, dt, noise_on, 0, pars);
    Z_temp(i,:,1:length(t)) = Z;
end    

Z_12(:,:) = Z_temp(1,:,1:length(0:dt:12));
Z_18(:,:) = Z_temp(2,:,1:length(0:dt:18));
Z_24(:,:) = Z_temp(3,:,1:length(0:dt:24));
Z_30(:,:) = Z_temp(4,:,1:length(0:dt:30));
Z_36(:,:) = Z_temp(5,:,1:length(0:dt:36));
Z_42(:,:) = Z_temp(6,:,1:length(0:dt:42));
Z_48(:,:) = Z_temp(7,:,1:length(0:dt:48));
Z_72(:,:) = Z_temp(8,:,1:length(0:dt:72));



dt = 20/3600;    % dt = 20s, tf = 12 hours
    
figure;
subplot(2,8,[1,2])
t = 0:dt:12;
plot(t, (Z_12(2,:)), t, (Z_12(3,:)), t, (Z_12(4,:)),t, (Z_12(6,:)),'Linewidth', line_thickness)
%xlabel('$Time~[hours]$','Interpreter','latex')
%ylabel('$Population$','Interpreter','latex')
xlim([0,12]);
xticks([0  6 12])
ylim([1e3, 1e9]);
yticks([1e3,1e5, 1e7, 1e9])
set(gca,'yscale','log')
yyaxis right
hold on
plot(t, Z_12(7,:)*conversion_factor,'Linewidth', line_thickness);
ylim([1e0, 1e4]);
yticks([1e0,1e2, 1e4])
title('$\mathrm{T_{max}=12}$','Interpreter','latex')
set(gca,'FontSize',20);
fig.PaperUnits = 'inches';
%pbaspect([2.5 1.5 1])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')


%figure;
subplot(2,8,[9,11])
t = 0:dt:18;
plot(t, (Z_18(2,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_18(3,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_18(4,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_18(6,:)),  'Linewidth', line_thickness)
hold on
%xlabel('$Time~[hours]$','Interpreter','latex')
%ylabel('$Population$','Interpreter','latex')

title('$\mathrm{T_{max}=18}$','Interpreter','latex')
xlim([0,18]);
xticks([0  6 12 18])
ylim([1e3, 1e9]);
yticks([1e3,1e5, 1e7, 1e9])
%legend('S','E','L','V')
set(gca,'yscale','log')
yyaxis right
hold on
plot(t, Z_18(7,:)*conversion_factor,'Linewidth', line_thickness);
ylim([1e0, 1e4]);
yticks([1e0,1e2, 1e4])
set(gca,'FontSize',20);
fig.PaperUnits = 'inches';
%pbaspect([2.5 1.5 1])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')

figure;
subplot(2,8,[1,4]);
t = 0:dt:24;
plot(t, (Z_24(2,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_24(3,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_24(4,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_24(6,:)),  'Linewidth', line_thickness)
hold on
%xlabel('$Time~[hours]$','Interpreter','latex')
%ylabel('$\mathrm{Population~[ml^{-1}]}$','Interpreter','latex')
set(gca,'yscale','log')

title('$\mathrm{T_{max}=24}$','Interpreter','latex')
xlim([0,24]);
xticks([0  6 12 18 24])
ylim([1e3, 1e9]);
yticks([1e3 1e5 1e7])
yyaxis right
hold on
plot(t, Z_24(7,:)*conversion_factor,'Linewidth', line_thickness);
ylim([1e0, 1e4]);
yticks([1e0,1e2, 1e4])
ylabel('$\mathrm{Arbitrium~[nM]}$','Interpreter','latex')

%legend('S','E','L','V')
set(gca,'FontSize',20);
fig.PaperUnits = 'inches';
%pbaspect([2.5 1.5 1])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')

figure;
subplot(2,8,[1,6]);
t = 0:dt:36;
plot(t, (Z_36(2,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_36(3,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_36(4,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_36(6,:)),  'Linewidth', line_thickness)
hold on
%xlabel('$Time~[hours]$','Interpreter','latex')
%ylabel('$Population$','Interpreter','latex')
xlim([0,36]);
xticks([0  12 24 36])
ylim([1e3, 1e9]);
yticks([1e3 1e5 1e7 ])
%legend('S','E','L','V')
set(gca,'yscale','log')
yyaxis right
hold on
plot(t, Z_36(7,:)*conversion_factor,'Linewidth', line_thickness);
ylim([1e0, 1e4]);
yticks([1e0,1e2, 1e4])
title('$\mathrm{T_{max}=36}$','Interpreter','latex')
set(gca,'FontSize',20);
fig.PaperUnits = 'inches';
%pbaspect([2.5 1.5 1])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')


%figure;
subplot(2,8,[9,16]);
t = 0:dt:48;
plot(t, (Z_48(2,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_48(3,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_48(4,:)),  'Linewidth', line_thickness)
hold on
plot(t, (Z_48(6,:)),  'Linewidth', line_thickness)
hold on
xlabel('Time [hours]','Interpreter','latex')
%ylabel('$Population$','Interpreter','latex')
xlim([0,48]);
xticks([0  12 24 36 48])
ylim([1e3, 1e9]);
yticks([1e3 1e5 1e7])
%legend('S','E','L','V')
set(gca,'yscale','log')
yyaxis right
hold on
plot(t, Z_48(7,:)*conversion_factor,'Linewidth', line_thickness);
ylim([1e0, 1e4]);
yticks([1e0,1e2, 1e4])
title('$\mathrm{T_{max}=48}$','Interpreter','latex')
set(gca,'FontSize',20);
fig.PaperUnits = 'inches';
%pbaspect([2.5 1.5 1])
legend('S','E','L','V','A','Interpreter', 'latex')
legend boxoff
set(gca,'yscale','log')
set(gca,'TickLabelInterpreter','latex')

%% E(tf) and L(tf)
state_final = zeros(2,length(tf_vector));
frac_end = state_final;
for tf_idx = 1:length(tf_vector)
    dt = 20/3600;    % dt = 20s, tf = 12 hours
    t = 0:dt:tf_vector(tf_idx);
    Z = forward_euler(Z0, optimal_params(:,tf_idx), t, dt, noise_on, 0, pars);
    state_final(:,tf_idx) = Z(3:4,end);
    frac_end(1,tf_idx) = state_final(1,tf_idx)/(state_final(1,tf_idx)+state_final(2,tf_idx));
    frac_end(2,tf_idx) = 1 - frac_end(1,tf_idx);
end    
%%  E(t) and L(t) for fixed vs optimal
mixed_prob = 0.5;
% tf_vals = [4,10];
tf_vals = 1:length(tf_vector);
state_final_optimal = zeros(2,length(tf_vector));
state_final_lysis = zeros(2,length(tf_vector));
state_final_lysogeny = zeros(2,length(tf_vector));
state_final_half = zeros(2,length(tf_vector));

for tf_idx = 1:length(tf_vals)
    
    dt = 20/3600;    % dt = 20s, tf = 12 hours
    t = 0:dt:tf_vector(tf_vals(tf_idx));
    Z = forward_euler(Z0, optimal_params(:,tf_vals(tf_idx)), t, dt, noise_on, 0, pars);
        
    Z_fixed_lysis = forward_euler(Z0, optimal_params(:,tf_vals(tf_idx)), t, dt, noise_on, 0.01, pars);
    Z_fixed_mixed = forward_euler(Z0, optimal_params(:,tf_vals(tf_idx)), t, dt, noise_on, 0.5, pars);
    Z_fixed_lysogeny = forward_euler(Z0, optimal_params(:,tf_vals(tf_idx)), t, dt, noise_on, 1, pars);
    
    state_final_optimal(:,tf_idx) = Z(3:4,end);
    state_final_lysis(:,tf_idx) = Z_fixed_lysis(3:4,end);
    state_final_lysogeny(:,tf_idx) = Z_fixed_lysogeny(3:4,end);
    state_final_half(:,tf_idx) = Z_fixed_mixed(3:4,end);
end

%% Fig. 5(a)
Z_end = zeros(12,3);
t = 0:dt:48;
idx = 1;
for prob = 0:0.1:1
    idx = idx+1;
    prob_new = prob;
    if prob == 0
        prob_new = 0.01;
    end    
    Z_temp = forward_euler(Z0, optimal_params(:,7), t, dt, noise_on, prob_new, pars);
    Z_end(idx,1:2) = Z_temp(3:4,end);
    Z_end(idx,3) = Z_end(idx,1)+Z_end(idx,2);
end

Z_end(1,1:2) = Z_48(3:4,end);
Z_end(1,3) = Z_end(1,1)+Z_end(1,2);
Z_end_log = log(Z_end);
figure;
H = bar((Z_end(:,1:2)),'stacked');
ylabel('$\mathrm{Population~[ml^{-1}]}$','Interpreter','latex')
xlabel('Lysogeny Probability, P','Interpreter','latex')
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gca,'FontSize',30);
set(gca,'yscale','log')
xticks(1:12)
yticks([1,1e2,1e4,1e6])

xticklabels({'Optimal', 'P=0','P=0.1','P=0.2','P=0.3','P=0.4','P=0.5','P=0.6','P=0.7','P=0.8','P=0.9','P=1'})
xticklabels({'Optimal', '0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
legend('Exposed Cells','Lysogens','Interpreter','latex')
H(1).FaceColor = [0.3,0.3,0.3];
H(2).FaceColor = [0.7,0.7,0.7];
axis tight
ylim([1,1e6])
yticks([1,1e2,1e4,1e6])
set(gca,'TickLabelInterpreter','latex')
pbaspect([2.5 1.5 1])
legend boxoff

%% Figure 5(b) birth states comparison with fixed probability

figure;
plot(tf_vector, (state_final_optimal(1,:)+state_final_optimal(2,:)), 'color','k', 'Marker',ms{1}, 'MarkerSize',9,'Linewidth', line_thickness,'LineStyle', ls{1})
hold on
plot(tf_vector, (state_final_lysis(1,:)+state_final_lysis(2,:)), 'color','k', 'Marker',ms{2},'MarkerSize',9,'Linewidth', line_thickness,'LineStyle', ls{2})
hold on
plot(tf_vector, (state_final_half(1,:)+state_final_half(2,:)), 'color','k', 'Marker',ms{3}, 'MarkerSize',9,'Linewidth', line_thickness,'LineStyle', ls{3})
hold on    
plot(tf_vector, (state_final_lysogeny(1,:)+state_final_lysogeny(2,:)),'color','k', 'Marker',ms{4}, 'MarkerSize',9,'Linewidth', line_thickness,'LineStyle', ls{4})
hold on 
xlabel('Time [hours]','Interpreter','latex')
ylabel('$\mathrm{Population~[ml^{-1}]}$','Interpreter','latex')
ylim([1e2,1e7])
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8])
xlim([12,48])
xticks([12 24 36 48])
fig.PaperUnits = 'inches';
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gca,'FontSize',30);
fig.PaperUnits = 'inches';
pbaspect([2.5 1.5 1])
legend('E+L (Optimal)','E+L (Lysis)', 'E+L (Mixed)','E+L (Lysogeny)','Interpreter','latex')
set(gca,'yscale','log')
legend boxoff
set(gca,'TickLabelInterpreter','latex')




%% subsequent code for Fig. 6
%% create matrix of optimal params for all J and r0
optimal_params_different_J = zeros(4,2,8);
for i = 0:3
    load (['optimal_params_r0_40_j_' num2str(i) '.mat']);
    optimal_params_different_J(i+1,:,:) = optimal_params;
end    
%save('optimal_params_different_J.mat','optimal_params_different_J');

optimal_params_different_r0 = zeros(4,2,8);
for i = 1:4
    load (['optimal_params_r0_' num2str(20+20*i) '_j_0.mat']);
    optimal_params_different_r0(i,:,:) = optimal_params;
end 
%save('optimal_params_different_r0.mat','optimal_params_different_r0');


phi = 3.4*10^(-10);
d_V = 0.05; 
E_0 = phi*S0*V0/(phi*S0+d_V);

%% heat map r0

%% Fig. 6(a) switch point of r0
figure;
matrix_vals = zeros(4,8);
for i = 1:4
    for j = 1:8
        matrix_vals(5-i,j) = optimal_params_different_r0(i,2,j)*conversion_factor*1e12;
    end
end    

cdata = matrix_vals(:,1:7);
yvalues = {'R_0 = 40','R_0 = 60','R_0 = 80','R_0 = 100'};
xvalues = {'12','18','24','30','36','42','48'};
%[gX gY]= meshgrid(12:6:48,40:20:100);

h = heatmap(xvalues,yvalues,cdata,'CellLabelColor','none','ColorLimits',[0 2000]);
set(gca,'FontSize',20);

h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
a2 = axes('Position', h.Position);               %new axis on top
a2.Color = 'none';                               %new axis transparent
a2.YTick = 1:size(h.ColorData,1);                %set y ticks to number of rows
a2.XTick = 1:size(h.ColorData,2);                %Remove xtick
ylim(a2, [0.5, size(h.ColorData,1)+.5])          %center the tick marks

a2.YTickLabel = {'$R_0=40$','$R_0=60$','$R_0=80$','$R_0=100$'}; 

xlim(a2, [0.5, size(h.ColorData,2)+.5])          %center the tick marks
a2.XTickLabel = {'$12$','$18$','$24$','$30$','$36$','$42$','$48$'}; 
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',20);

%% compute end state r0

L_matrix_diff_r0 = zeros(4,8);
E_matrix_diff_r0 = zeros(4,8);


for i = 1:4
    Z0(1) = 20 + 20*i;
for j = 1:8
    optimal_param = [optimal_params_different_r0(i,1,j),optimal_params_different_r0(i,2,j)];
    t = 0:dt:tf_vector(j);
    Z = forward_euler(Z0, optimal_param, t, dt, noise_on, 0, pars);
    E_matrix_diff_r0(i,j) = Z(3,end);
    L_matrix_diff_r0(i,j) = Z(4,end);
end
end
Z0(1) = 40;

fitness_matrix_r0 = L_matrix_diff_r0 + E_matrix_diff_r0;
L_frac_matrix_r0 = L_matrix_diff_r0./(L_matrix_diff_r0 + E_matrix_diff_r0);

L_frac_matrix_r0_saved = L_frac_matrix_r0;
fitness_matrix_r0_saved = fitness_matrix_r0;
%save('L_frac_matrix_r0_saved.mat','L_frac_matrix_r0_saved')   %(only do once)
%save('fitness_matrix_r0_saved.mat','fitness_matrix_r0_saved') %(only do once)
%% Fig. 6(b) growth rate r0
load('L_frac_matrix_r0_saved.mat')
load('fitness_matrix_r0_saved.mat')

figure;
matrix_vals = zeros(4,8);
for i = 1:4
    for j = 1:8
        matrix_vals(5-i,j) = L_frac_matrix_r0_saved(i,j);
    end
end    

cdata = matrix_vals(:,1:7);
yvalues = {'J = 0','J = 1','J = 2','J = 3'};
xvalues = {'12','18','24','30','36','42','48'};
%[gX gY]= meshgrid(12:6:48,40:20:100);

h = heatmap(xvalues,yvalues,cdata,'CellLabelColor','none','ColorLimits',[0.7 1]);
set(gca,'FontSize',20);

h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
a2 = axes('Position', h.Position);               %new axis on top
a2.Color = 'none';                               %new axis transparent
a2.YTick = 1:size(h.ColorData,1);                %set y ticks to number of rows
a2.XTick = 1:size(h.ColorData,2);                %Remove xtick
ylim(a2, [0.5, size(h.ColorData,1)+.5])          %center the tick marks

a2.YTickLabel = {'$R_0=40$','$R_0=60$','$R_0=80$','$R_0=100$'}; 

xlim(a2, [0.5, size(h.ColorData,2)+.5])          %center the tick marks
a2.XTickLabel = {'$12$','$18$','$24$','$30$','$36$','$42$','$48$'}; 


set(gca,'TickLabelInterpreter','latex')

set(gca,'FontSize',20);

%% Fig. 6(c) frac_L r0
figure;
matrix_vals = zeros(4,8);
growth_rate = zeros(4,8);

for i = 1:4
    for j = 1:8
        growth_rate(i,j) = (1/tf_vector(j))*log(fitness_matrix_r0_saved(i,j)/E_0);
        matrix_vals(5-i,j) = growth_rate(i,j);
    end
end    

cdata = matrix_vals(:,1:7);
yvalues = {'J = 0','J = 1','J = 2','J = 3'};
xvalues = {'12','18','24','30','36','42','48'};
%[gX gY]= meshgrid(12:6:48,40:20:100);

h = heatmap(xvalues,yvalues,cdata,'CellLabelColor','none','ColorLimits',[0 0.7]);
set(gca,'FontSize',20);

h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
a2 = axes('Position', h.Position);               %new axis on top
a2.Color = 'none';                               %new axis transparent
a2.YTick = 1:size(h.ColorData,1);                %set y ticks to number of rows
a2.XTick = 1:size(h.ColorData,2);                %Remove xtick
ylim(a2, [0.5, size(h.ColorData,1)+.5])          %center the tick marks

a2.YTickLabel = {'$R_0=40$','$R_0=60$','$R_0=80$','$R_0=100$'}; 

xlim(a2, [0.5, size(h.ColorData,2)+.5])          %center the tick marks
a2.XTickLabel = {'$12$','$18$','$24$','$30$','$36$','$42$','$48$'}; 
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',20);

%% Fig. 6(d) optimal params J
figure;
matrix_vals = zeros(4,8);
for i = 1:4
    for j = 1:8
        matrix_vals(5-i,j) = optimal_params_different_J(i,2,j)*conversion_factor*1e12;
    end
end    

cdata = matrix_vals(:,1:7);
yvalues = {'J = 0','J = 1','J = 2','J = 3'};
xvalues = {'12','18','24','30','36','42','48'};
%[gX gY]= meshgrid(12:6:48,40:20:100);

h = heatmap(xvalues,yvalues,cdata,'CellLabelColor','none','ColorLimits',[0 2000]);
set(gca,'FontSize',20);

h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
a2 = axes('Position', h.Position);               %new axis on top
a2.Color = 'none';                               %new axis transparent
a2.YTick = 1:size(h.ColorData,1);                %set y ticks to number of rows
a2.XTick = 1:size(h.ColorData,2);                %Remove xtick
ylim(a2, [0.5, size(h.ColorData,1)+.5])          %center the tick marks

a2.YTickLabel = {'$J = 0$','$J = 1$','$J = 2$','$J = 3$'}; 

xlim(a2, [0.5, size(h.ColorData,2)+.5])          %center the tick marks
a2.XTickLabel = {'$12$','$18$','$24$','$30$','$36$','$42$','$48$'}; 
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',20);

%% compute end state
L_matrix_diff_J = zeros(4,8);
E_matrix_diff_J = zeros(4,8);

for i = 1:4
    pars.J = i-1;
for j = 1:8
    optimal_param = [optimal_params_different_J(1,1,j),optimal_params_different_J(1,2,j)];
    t = 0:dt:tf_vector(j);
    Z = forward_euler(Z0, optimal_param, t, dt, noise_on, 0, pars);
    E_matrix_diff_J(i,j) = Z(3,end);
    L_matrix_diff_J(i,j) = Z(4,end);
end
end

fitness_matrix_J = L_matrix_diff_J + E_matrix_diff_J;
L_frac_matrix_J = L_matrix_diff_J./(L_matrix_diff_J + E_matrix_diff_J);
  
L_frac_matrix_J_saved = L_frac_matrix_J;
fitness_matrix_J_saved = fitness_matrix_J;
%save('L_frac_matrix_J_saved.mat','L_frac_matrix_J_saved')   %(only do once)
%save('fitness_matrix_J_saved.mat','fitness_matrix_J_saved') %(only do once)

%% Fig 6(e) frac_L J

load('L_frac_matrix_J_saved.mat')
load('fitness_matrix_J_saved.mat')

figure;
matrix_vals = zeros(4,8);
for i = 1:4
    for j = 1:8
        matrix_vals(5-i,j) = L_frac_matrix_J_saved(i,j);
    end
end    

cdata = matrix_vals(:,1:7);
yvalues = {'J = 0','J = 1','J = 2','J = 3'};
xvalues = {'12','18','24','30','36','42','48'};
%[gX gY]= meshgrid(12:6:48,40:20:100);

h = heatmap(xvalues,yvalues,cdata,'CellLabelColor','none','ColorLimits',[0.7 1]);
set(gca,'FontSize',20);

h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
a2 = axes('Position', h.Position);               %new axis on top
a2.Color = 'none';                               %new axis transparent
a2.YTick = 1:size(h.ColorData,1);                %set y ticks to number of rows
a2.XTick = 1:size(h.ColorData,2);                %Remove xtick
ylim(a2, [0.5, size(h.ColorData,1)+.5])          %center the tick marks

a2.YTickLabel = {'$J = 0$','$J = 1$','$J = 2$','$J = 3$'}; 
xlim(a2, [0.5, size(h.ColorData,2)+.5])          %center the tick marks
a2.XTickLabel = {'$12$','$18$','$24$','$30$','$36$','$42$','$48$'}; 

set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',20);

%% Fig. 6(f) growth rate J
figure;
matrix_vals = zeros(4,8);
growth_rate = zeros(4,8);

for i = 1:4
    for j = 1:8
        growth_rate(i,j) = (1/tf_vector(j))*log(fitness_matrix_J_saved(i,j)/E_0);
        matrix_vals(5-i,j) = growth_rate(i,j);
    end
end    

cdata = matrix_vals(:,1:7);
yvalues = {'J = 0','J = 1','J = 2','J = 3'};
xvalues = {'12','18','24','30','36','42','48'};
%[gX gY]= meshgrid(12:6:48,40:20:100);

h = heatmap(xvalues,yvalues,cdata,'CellLabelColor','none','ColorLimits',[0 0.7]);
set(gca,'FontSize',20);

h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
a2 = axes('Position', h.Position);               %new axis on top
a2.Color = 'none';                               %new axis transparent
a2.YTick = 1:size(h.ColorData,1);                %set y ticks to number of rows
a2.XTick = 1:size(h.ColorData,2);                %Remove xtick
ylim(a2, [0.5, size(h.ColorData,1)+.5])          %center the tick marks

a2.YTickLabel = {'$J = 0$','$J = 1$','$J = 2$','$J = 3$'}; 

xlim(a2, [0.5, size(h.ColorData,2)+.5])          %center the tick marks
a2.XTickLabel = {'$12$','$18$','$24$','$30$','$36$','$42$','$48$'}; 


set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',20);

%% figure S1 prob vs A for different production rates
load('optimal_params_r0_40_j_0_sweep.mat')
conversion_factor = 10^3*10^9/(6.022*10^23);    % ml, nano, mole

figure;
prod_rate_vector_log = 5:0.25:8;
prod_rate_vector = 10.^prod_rate_vector_log;
plot(prod_rate_vector, conversion_factor*10^12*optimal_params(1:13,2,3), 'color','k', 'Marker',ms{1}, 'MarkerSize',9,'Linewidth', line_thickness,'LineStyle', ls{1})

xlabel('Generation~rate, $k_{L, I}$ ~$\mathrm{[molecules/cell~h^{-1}]}$','Interpreter','latex')
ylabel('Switching~concentration [nM]','Interpreter','latex')
xlim([10^5,10^8]);
set(gca,'TickLabelInterpreter','latex')
fig.PaperUnits = 'inches';
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gca,'FontSize',25);
set(gca,'xscale','log')
set(gca,'yscale','log')

fig.PaperUnits = 'inches';
pbaspect([2.5 1.5 1])

xvals = prod_rate_vector;
yvals = conversion_factor*10^12*optimal_params(1:13,2,3);

Sxx = sum(xvals.^2) - sum(xvals)^2/length(xvals);
Syy = sum(yvals.^2) - sum(yvals)^2/length(yvals);
Sxy = dot(xvals, yvals) - sum(xvals)*sum(yvals)/length(yvals);

R = Sxy/(sqrt(Sxx*Syy));
