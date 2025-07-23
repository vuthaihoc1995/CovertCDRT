clc;
clear variables;
close all;
%% Main parameter
PdBm = 0:2.5:30;
P = db2pow(PdBm);

sigma = 1;

alpha_j = 0.7;
alpha_i = 1-alpha_j;
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


    %% fixed alpha_j = 0.7;

   snr_RUi = @(x,y) P(idx)*x.*y/(P(idx)*(1-x).*y + sigma);
   snr_Uii = @(x) P(idx)*(1-x).*mean(abs(h_SUi).^2)/sigma;
   % fixed alpha_i = 0.7;
   C_ui1_fixed (idx) =1/2*log2( 1 + snr_Uii(alpha_j)); 
   R_uij_fixed (idx) =1/2*log2( 1 + snr_RUi(alpha_j,mean(abs(h_SUi).^2))); 
   R_r_fixed (idx) =1/2*log2( 1 + snr_RUi(alpha_j,mean(abs(h_SR).^2))); 


  %% Optimal solution ;
  epsilon = 0.99; Rj = 0.5;
  %%%% Find root of epsilon[1];
  fun_alpha = @(x) 1 - x.^(-x./(x-1)) + x.^(-1./(x-1)) - epsilon;
  x0 = (1-0.5).*rand(1,1) + 0.5;
  alphaj_t = fzero(fun_alpha,x0);

  %%% Feasible solution
  zeta = max(    mean(abs(h_SR).^2) , mean(abs(h_SUi).^2));
  % cond = ( P(idx)*zeta +sigma ) /( P(idx)*zeta)*(1- 1/2^(2*Rj))
  
  cond = 1  - (  (P(idx)*zeta +sigma )/2^(2*Rj)   - 1 )   /(P(idx)*zeta);

  xi = max(cond,1/2);

  % Final solution;
  alpha_opt(idx) = min(1,min(alphaj_t,xi));

   C_ui1_opt (idx) = 1/2*log2( 1 + snr_Uii(alpha_opt(idx))); 
    R_uij_opt (idx) =1/2*log2( 1 + snr_RUi(alpha_opt(idx),mean(abs(h_SUi).^2))); 
    R_r_opt (idx) =1/2*log2( 1 + snr_RUi(alpha_opt(idx),mean(abs(h_SR).^2))); 
        
end

Fix1 = plot(PdBm,C_ui1_fixed,'k-s','MarkerSize',8,'LineWidth',1.5); hold on;
Fix2 = plot(PdBm,R_uij_fixed,'y-x','MarkerSize',8,'LineWidth',1.5); hold on;
Fix3 = plot(PdBm,R_r_fixed,'b-v','MarkerSize',8,'LineWidth',1.5); hold on;
Fix4 =yline(fun_alpha(alpha_j)+ epsilon,'k-x','LineWidth',1.5); hold on;

req1 = yline(Rj,'m--','LineWidth',1.5); hold on;
req2 = yline(epsilon,'g-.','LineWidth',1.5); hold on;


Opt1 = plot(PdBm,C_ui1_opt,'k-o','MarkerSize',10,'LineWidth',1.5); hold on;
Opt2 = plot(PdBm,R_uij_opt,'b--+','MarkerSize',8,'LineWidth',1.5); hold on;
Opt3 = plot(PdBm,R_r_opt,'r-^','MarkerSize',9,'LineWidth',1.5); hold on;
% Opt4 =yline(fun_alpha(alpha_opt)+ epsilon,'r-*','LineWidth',1.5); hold on;


hold on;
colormap(gca) % cool

lgd=legend([Fix1(1), Fix2(1), Fix3(1),Fix4(1),...
    Opt1(1), Opt2(1), Opt3(1),...
    req1(1),req2(1)],...
    '$\alpha_j = 0.7$: $C_{i}[1]$','$\alpha_j = 0.7$: $R_{\rm R}[1]$','$\alpha_j = 0.7$: $R_{{\rm U}_i}[1]$','$\alpha_j = 0.7 $: $\epsilon_{\rm E}[1]$',...
     'Eq. (17): $C_{i}[1]$','Eq. (17): $R_{\rm R}[1]$','Eq. (17): $R_{{\rm U}_i}[1]$ ',...
    'Rate: $R_j=0.5$',...
    'DEP: $\epsilon_{\rm th}=0.9$'...
    ,'FontSize',11,'location','se','Interpreter','latex');
xlabel('Transmit Power, {\it P} [dBm]','Fontsize',16) 
ylabel('Covert Rate [bps/Hz]');
axis([min(PdBm) max(PdBm) 0 3]);
lgd.NumColumns = 1;
lgd.FontSize = 14;
set(gca,'fontsize',14);

% saveas(gcf,'Figure3a.fig');

