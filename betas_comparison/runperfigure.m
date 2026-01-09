mvecs = 110:100:1010;

load("results_betas1.mat")
colors = [ 0.36 0.14 0.32; 
           1.0 0.5 0.0;
           0.0 0.8 0.0;
           1 0 1;
           0.0 0.3 1.0;
           ];

figure(1)

subplot(2,3,1)
hold on
errorbar(mvecs-10, beta18./trG, SD_beta18./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
errorbar(mvecs-5, beta14./trG, SD_beta14./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs, beta12./trG, SD_beta12./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs+5, beta34./trG, SD_beta34./trG, '-', 'Color', colors(4,:), 'LineWidth', 4);
errorbar(mvecs+10, beta78./trG, SD_beta78./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Alg','fontsize',15)
%legend('error (A)', 'bound (A)', 'error (B)','bound (B)', 'error (C)', 'bound (C)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


load("results_betas2.mat")

figure(1)

subplot(2,3,2)
hold on
errorbar(mvecs-10, beta18./trG, SD_beta18./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
errorbar(mvecs-5, beta14./trG, SD_beta14./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs, beta12./trG, SD_beta12./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs+5, beta34./trG, SD_beta34./trG, '-', 'Color', colors(4,:), 'LineWidth', 4);
errorbar(mvecs+10, beta78./trG, SD_beta78./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Geom','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


load("results_betas3.mat")


figure(1)

subplot(2,3,3)
hold on
errorbar(mvecs-10, beta18./trG, SD_beta18./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
errorbar(mvecs-5, beta14./trG, SD_beta14./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs, beta12./trG, SD_beta12./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs+5, beta34./trG, SD_beta34./trG, '-', 'Color', colors(4,:), 'LineWidth', 4);
errorbar(mvecs+10, beta78./trG, SD_beta78./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Gaps','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 



load("results_betas4.mat")

figure(1)

subplot(2,3,4)
hold on
errorbar(mvecs-10, beta18./trG, SD_beta18./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
errorbar(mvecs-5, beta14./trG, SD_beta14./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs, beta12./trG, SD_beta12./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs+5, beta34./trG, SD_beta34./trG, '-', 'Color', colors(4,:), 'LineWidth', 4);
errorbar(mvecs+10, beta78./trG, SD_beta78./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example RBF','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 



load("results_betas5.mat")

figure(1)

subplot(2,3,5)
hold on
errorbar(mvecs-10, beta18./trG, SD_beta18./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
errorbar(mvecs-5, beta14./trG, SD_beta14./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs, beta12./trG, SD_beta12./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs+5, beta34./trG, SD_beta34./trG, '-', 'Color', colors(4,:), 'LineWidth', 4);
errorbar(mvecs+10, beta78./trG, SD_beta78./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(1/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


load("results_betas6.mat")

figure(1)

subplot(2,3,6)
hold on
h1 = errorbar(mvecs-10, beta18./trG, SD_beta18./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
h2 = errorbar(mvecs-5, beta14./trG, SD_beta14./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
h3 = errorbar(mvecs, beta12./trG, SD_beta12./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
h4 = errorbar(mvecs+5, beta34./trG, SD_beta34./trG, '-', 'Color', colors(4,:), 'LineWidth', 4);
h5 = errorbar(mvecs+10, beta78./trG, SD_beta78./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(3/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 

ax_leg = axes('Position',[0.1 0.01 0.8 0.05],'Visible','off'); 
lgd = legend(ax_leg, [h1 h2 h3 h4 h5], {'\beta=1/8', '\beta=1/4', '\beta=1/2', '\beta=3/4','\beta=7/8'}, ...
    'Orientation','horizontal','FontSize',30,'NumColumns',5);
lgd.Box = 'off';         % rimuove il bordo della legenda

