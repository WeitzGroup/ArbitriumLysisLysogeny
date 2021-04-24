prob12 = zeros(1,length(t));
prob18 = zeros(1,length(t));
prob24 = zeros(1,length(t));
prob48 = zeros(1,length(t));
prob72 = zeros(1,length(t));
for i =1:length(t)
    prob12(i) = probability(Z_12(7,i),optimal_params_save(:,1), 10^12,1, 0);
end    

for i =1:length(t)
    prob18(i) = probability(Z_12(7,i),optimal_params_save(:,1), 10^12,1, 0);
end
for i =1:length(t)
    prob24(i) = probability(Z_12(7,i),optimal_params_save(:,1), 10^12,1, 0);
end
for i =1:length(t)
    prob48(i) = probability(Z_12(7,i),optimal_params_save(:,1), 10^12,1, 0);
end
for i =1:length(t)
    prob72(i) = probability(Z_12(7,i),optimal_params_save(:,1), 10^12,1, 0);
end
