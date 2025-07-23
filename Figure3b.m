clc;
clear variables;
close all;
%% Main parameter
PdBm = 0:2.5:30;
P = db2pow(PdBm);

sigma = 1;

beta_s = 0.7;


sim_times = 1;

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
epsilon = 0.99; Rj = 0.5;
for idx = 1:length(PdBm)
    
    %% Fixed power 
    snr_Uj  = @(x) (1-x)*P(idx).*mean(abs(h_RUj).^2)/ sigma;
    snr_Ui = @(x) x*P(idx).*mean(abs(h_SUi).^2)/ sigma;

    C_ui_fixed(idx) = 1/2*log2( 1 + snr_Ui(beta_s));   
    R_uj_fixed(idx) = 1/2*log2( 1 + snr_Uj(beta_s));   

    Delta= @(x) lambda_RE/lambda_SE*(1./x -1);
    fun_beta = @(x) 1 - 1./(Delta(x) - 1 ).*(Delta(x).^(  1/(  1-Delta(x)  ) )  - Delta(x).^(   Delta(x)/(1-Delta(x))  )   ) - epsilon;

    %%
   %%%% Find root of epsilon[2];
  x0 = 0.9;
  beta_t = fzero(fun_beta,x0);

  %%% Feasible solution
  Phi = min(1,1- (2^(2*Rj)  - 1 )*sigma/(P(idx)*mean(abs(h_RUj).^2) )  );


  % Final solution;
  beta_opt = max(1/2,min(beta_t,Phi));
  C_ui_opt(idx) = 1/2*log2( 1 + snr_Ui(beta_opt));   
R_uj_opt(idx) = 1/2*log2( 1 + snr_Uj(beta_opt));   
        
end

Fix1 = plot(PdBm,C_ui_fixed,'k-s','MarkerSize',8,'LineWidth',1.5); hold on;
Fix2 = plot(PdBm,R_uj_fixed,'b-x','MarkerSize',8,'LineWidth',1.5); hold on;
Fix3 = yline(fun_beta(beta_s)+ epsilon,'k-.','LineWidth',1.5); hold on;

req1 = yline(Rj,'m-','LineWidth',1.5); hold on;
req2 = yline(epsilon,'g--','LineWidth',1.5); hold on;

Opt1 = plot(PdBm,C_ui_opt,'k-o','MarkerSize',10,'LineWidth',1.5); hold on;
Opt2 = plot(PdBm,R_uj_opt,'r-v','MarkerSize',8,'LineWidth',1.5); hold on;

hold on;
colormap(gca) % cool

lgd=legend([Fix1(1), Fix2(1), Fix3(1),...
    Opt1(1), Opt2(1),...
    req1(1),req2(1)],...
    '$\beta_s = 0.7$: $C_{{\rm U}_i}[2]$','$\beta_s = 0.7$: $R_{{\rm U}_j}[2]$','$\beta_s = 0.7$: $\epsilon_{\rm E}[2]$',...
     'Eq. (19): $C_{{\rm U}_i}[2]$','Eq. (19): $R_{{\rm U}_j}[2]$ ',...
    'Rate: $R_j=0.5$',...
    'DEP: $\epsilon_{\rm th}=0.99$'...
    ,'FontSize',11,'location','se','Interpreter','latex');
xlabel('Transmit Power, {\it P} [dBm]','Fontsize',16) 
ylabel('Covert Rate [bps/Hz]');
axis([min(PdBm) max(PdBm) 0 3]);
lgd.NumColumns = 1;
lgd.FontSize = 14;
set(gca,'fontsize',14);

% saveas(gcf,'Figure3b.fig');
