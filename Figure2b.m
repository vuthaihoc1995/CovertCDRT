clc;
clear variables;
close all;
%% Main parameter
P = db2pow(20);
sigma = 1;

tau_dBm = 0:2.5:30;
tau = db2pow(tau_dBm);

tau_dBm2 = 0:0.5:30;
tau2 = db2pow(tau_dBm2);

sim_times = 1e5;
%% Channel generation
lambda_SR = 1;
lambda_SUi = 0.5;

lambda_RUj = 0.5;
lambda_SE = 1;
lambda_RE = 0.5;

h_SR  = sqrt(lambda_SR/2)*(randn(1,sim_times) + 1i*randn(1,sim_times));
h_SUi = sqrt(lambda_SUi/2)*(randn(1,sim_times) + 1i*randn(1,sim_times));

h_RUj = sqrt(lambda_RUj/2)*(randn(1,sim_times) + 1i*randn(1,sim_times));
h_SE  = sqrt(lambda_SE/2)*(randn(1,sim_times) + 1i*randn(1,sim_times));
h_RE  = sqrt(lambda_RE/2)*(randn(1,sim_times) + 1i*randn(1,sim_times));

del = [0.5, 0.6, 0.7];
for ss = 1:length(del)
%% parameter 
alpha_j = del(ss);
alpha_i = 1-alpha_j;
beta_s = del(ss);
beta_r = 1-beta_s;


%% Simulation
for idx = 1:length(tau_dBm)
    % Phase 1
    calH0_1 = P*alpha_j.*abs(h_SE).^2 + sigma;
    calH1_1 = P*alpha_i.*abs(h_SE).^2 + P*alpha_j.*abs(h_SE).^2 + sigma;
    
    DEP_1_sim(ss,idx) = sum(calH0_1 > tau(idx))/sim_times + sum(calH1_1 < tau(idx))/sim_times;



    % Phase 2
    calH0_2 = beta_r*P.*abs(h_RE).^2 + sigma;
    calH1_2 = beta_s*P.*abs(h_SE).^2 + beta_r*P.*abs(h_RE).^2 + sigma;
    
    DEP_2_sim(ss,idx) = sum(calH0_2 > tau(idx))/sim_times + sum(calH1_2 < tau(idx))/sim_times;
    
end

%% Analytical
for ii = 1:length(tau2)
    % Phase 1
    if tau2(ii) - sigma < 0
        DEP_1_ana(ss,ii) = 1;
    else
        DEP_1_ana(ss,ii) = 1 - exp(-(tau2(ii)-sigma)/(P*lambda_SE)) + exp(-(tau2(ii)-sigma)/(P*alpha_j*lambda_SE));
    end


    % Phase 2
    if tau2(ii) - sigma < 0
        DEP_2_ana(ss,ii) = 1;
    else
        DEP_2_ana(ss,ii) = 1 - lambda_SE*beta_s/(lambda_RE*beta_r - lambda_SE*beta_s) * (exp(-(tau2(ii)-sigma)/(lambda_RE*beta_r*P)) - exp(-(tau2(ii)-sigma)/(lambda_SE*beta_s*P)));
    end
    
end
%% Optimal
tau_1_opt(ss) = pow2db(-P*lambda_SE*alpha_j/(1-alpha_j)*log(alpha_j) + sigma);
DEP_1_opt(ss) = 1 - alpha_j^(alpha_j/(1-alpha_j)) + alpha_j^(1/(1-alpha_j));
tau_2_opt(ss) = pow2db(-lambda_SE*beta_s*lambda_RE*beta_r*P/(lambda_SE*beta_s - lambda_RE*beta_r)*log(lambda_RE*beta_r/(lambda_SE*beta_s)) + sigma);
DEP_2_opt(ss) = 1 - lambda_SE*beta_s/(lambda_RE*beta_r - lambda_SE*beta_s) * ((lambda_RE*beta_r/(lambda_SE*beta_s))^(lambda_SE*beta_s/(lambda_SE*beta_s-lambda_RE*beta_r)) - (lambda_RE*beta_r/(lambda_SE*beta_s))^(lambda_RE*beta_r/(lambda_SE*beta_s-lambda_RE*beta_r)) );
end
%% Plot
blue1 = [0.00,0.45,0.74];  pink1 = [1.00,0.07,0.65];
green1 = [0.47,0.67,0.19]; orrange = [0.85,0.33,0.10];
% 
sim1 = plot(tau_dBm,DEP_1_sim,'ro','MarkerSize',10,'LineWidth',1.5,'Color',blue1); hold on;
ana1 = plot(tau_dBm2,DEP_1_ana,'k-','MarkerSize',10,'LineWidth',1.5); hold on;
opt1 = plot(tau_1_opt,DEP_1_opt,'mp','MarkerSize',10,'LineWidth',1.5); hold on;

sim2 = plot(tau_dBm,DEP_2_sim,'bs','MarkerSize',10,'LineWidth',1.5,'Color',orrange); hold on;
plot(tau_dBm2,DEP_2_ana,'k-','MarkerSize',10,'LineWidth',1.5); hold on;
plot(tau_2_opt,DEP_2_opt,'mp','MarkerSize',10,'LineWidth',1.5); hold on;

hold on;
colormap(gca) % cool

lgd=legend([sim1(1), sim2(1), ana1(1),opt1(1)],...
    'Simulation - Phase 1','Simulation - Phase 2',...
    'Theory - Eqs. (8), (11)',...
    'Optimality - Eqs. (10), (13)','FontSize',11,'location','se','Interpreter','latex');
xlabel('Judgment Threshold, \tau [dBm]','Fontsize',16) 
ylabel('Detection Error Probability (DEP)');
lgd.NumColumns = 1;
lgd.FontSize = 14;
set(gca,'fontsize',14);
axis([min(tau_dBm) max(tau_dBm) 0 1]);
saveas(gcf,'Figure2b.fig');

