clc;
clear variables;
close all;

lambda_SR = 1;
lambda_SUi = 0.5;

lambda_RUj = 0.5;
lambda_SE = 1;
lambda_RE = 0.5;

alpha = 0.5:0.1:1;

func1 = @(x) 1 - x.^(-x./(x-1)) + x.^(-1./(x-1));
Delta = @(x) lambda_RE*(1-x)./(lambda_SE*x);

func2 = @(x) 1 - 1./(Delta(x)-1) .* (Delta(x).^(1./(1-Delta(x)))   - Delta(x).^(Delta(x)./(1-Delta(x))));

%% Plot
blue1 = [0.00,0.45,0.74];  pink1 = [1.00,0.07,0.65];
green1 = [0.47,0.67,0.19]; orrange = [0.85,0.33,0.10];
% 
sim1 = plot(alpha,func1(alpha),'r-o','MarkerSize',10,'LineWidth',1.5); hold on;
sim2 = plot(alpha,func2(alpha),'b-s','MarkerSize',10,'LineWidth',1.5); hold on;
hold on;
colormap(gca) % cool
lgd=legend([sim1(1), sim2(1)],...
    'Optimality - Eq. (10)',...
    'Optimality - Eq. (13)','FontSize',11,'location','se','Interpreter','latex');
xlabel('Power Allocation Coefficient, \alpha_j/\beta_s','Fontsize',16) 
ylabel('Detection Error Probability (DEP)');
lgd.NumColumns = 1;
lgd.FontSize = 14;
set(gca,'fontsize',14);
axis([min(alpha) max(alpha) 0 1]);
saveas(gcf,'Figure2c.fig');

