clc;
clear variables;
close all;
%% Main parameter
PdBm = 0:2.5:30;
P = db2pow(PdBm);
PdBm2 = 0:0.5:30;
P2 = db2pow(PdBm2);

sigma = 1;

tau_dBm = 5;
tau = db2pow(tau_dBm);


alpha_j = 0.7;
alpha_i = 1-alpha_j;
beta_s = 0.5;
beta_r = 1-beta_s;


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


%% Simulation
for idx = 1:length(PdBm)


    % Phase 1
    calH0_1 = P(idx)*alpha_j.*abs(h_SE).^2 + sigma;
    calH1_1 = P(idx)*alpha_i.*abs(h_SE).^2 + P(idx)*alpha_j.*abs(h_SE).^2 + sigma;
    
    DEP_1_sim(idx) = sum(calH0_1 > tau)/sim_times + sum(calH1_1 < tau)/sim_times;



    % Phase 2
    calH0_2 = beta_r*P(idx).*abs(h_RE).^2 + sigma;
    calH1_2 = beta_s*P(idx).*abs(h_SE).^2 + beta_r*P(idx).*abs(h_RE).^2 + sigma;
    
    DEP_2_sim(idx) = sum(calH0_2 > tau)/sim_times + sum(calH1_2 < tau)/sim_times;
    
end

%% Analytical
for ii = 1:length(P2)

    % Phase 1
    if tau - sigma < 0
        DEP_1_ana(ii) = 1;
    else
        DEP_1_ana(ii) = 1 - exp(-(tau-sigma)/(P2(ii)*lambda_SE)) + exp(-(tau-sigma)/(P2(ii)*alpha_j*lambda_SE));
    end


    % Phase 2
    if tau - sigma < 0
        DEP_2_ana(ii) = 1;
    else
        DEP_2_ana(ii) = 1 - lambda_SE*beta_s/(lambda_RE*beta_r - lambda_SE*beta_s) * (exp(-(tau-sigma)/(lambda_RE*beta_r*P2(ii))) - exp(-(tau-sigma)/(lambda_SE*beta_s*P2(ii))));
    end
    

%% Optimal
tau_1_opt(ii) = pow2db(-P2(ii)*lambda_SE*alpha_j/(1-alpha_j)*log(alpha_j) + sigma);
tau_2_opt(ii) = pow2db(-lambda_SE*beta_s*lambda_RE*beta_r*P2(ii)/(lambda_SE*beta_s - lambda_RE*beta_r)*log(lambda_RE*beta_r/(lambda_SE*beta_s)) + sigma);

end

DEP_1_opt = 1 - alpha_j^(alpha_j/(1-alpha_j)) + alpha_j^(1/(1-alpha_j));
DEP_2_opt = 1 - lambda_SE*beta_s/(lambda_RE*beta_r - lambda_SE*beta_s) * ((lambda_RE*beta_r/(lambda_SE*beta_s))^(lambda_SE*beta_s/(lambda_SE*beta_s-lambda_RE*beta_r)) - (lambda_RE*beta_r/(lambda_SE*beta_s))^(lambda_RE*beta_r/(lambda_SE*beta_s-lambda_RE*beta_r)) );

%% Plot
blue1 = [0.00,0.45,0.74];  pink1 = [1.00,0.07,0.65];
green1 = [0.47,0.67,0.19]; orrange = [0.85,0.33,0.10];
% 
sim1 = plot(PdBm,DEP_1_sim,'ro','MarkerSize',10,'LineWidth',1.5,'Color',blue1); hold on;
ana1 = plot(PdBm2,DEP_1_ana,'k-','MarkerSize',10,'LineWidth',1.5); hold on;
opt1 = yline(DEP_1_opt,'m--','LineWidth',1.5); hold on;

sim2 = plot(PdBm,DEP_2_sim,'bs','MarkerSize',10,'LineWidth',1.5,'Color',orrange); hold on;
plot(PdBm2,DEP_2_ana,'k-','MarkerSize',10,'LineWidth',1.5); hold on;
yline(DEP_2_opt,'m--','LineWidth',1.5); hold on;

yyaxis right
thre1 = plot(PdBm2,tau_1_opt,'b-.','MarkerSize',10,'LineWidth',1.5); hold on;
plot(PdBm2,tau_2_opt,'b-.','MarkerSize',10,'LineWidth',1.5); hold on;

hold on;
colormap(gca) % cool

lgd=legend([sim1(1), sim2(1), ana1(1),thre1(1),opt1(1)],...
    'Simulation - Phase 1','Simulation - Phase 2',...
    'Theory - Eqs. (8), (11)',...
    'Theory - Eqs. (9), (12)',...
    'Optimality - Eqs. (10), (13)','FontSize',11,'location','se','Interpreter','latex');
xlabel('Transmit Power, {\it P} [dBm]','Fontsize',16) 
yyaxis left
ylabel('Detection Error Probability (DEP)');
axis([min(PdBm) max(PdBm) 0 1]);

yyaxis right
ylabel('Judgment Threshold, \tau [dBm]','Fontsize',16)
lgd.NumColumns = 1;
lgd.FontSize = 14;
set(gca,'fontsize',14);
axis([min(PdBm) max(PdBm) 0 30]);
saveas(gcf,'Figure2a.fig');

