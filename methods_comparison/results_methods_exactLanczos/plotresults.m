mvecs = 150:100:1050;

load("results_comparison1.mat")
colors = [ 0.36 0.14 0.32; 
           1.0 0.5 0.0; 
           0.0 0.8 0.0;
           1.0 0.0 1.0;
           0.0 0.3 1.0];

figure(1)

subplot(2,3,1)
hold on
errorbar(mvecs, err_detective./trG, SD_detective./trG, '-', 'Color', 'm', 'LineWidth', 4);
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Alg','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


load("results_comparison2.mat")

figure(1)

subplot(2,3,2)
hold on
errorbar(mvecs, err_detective./trG, SD_detective./trG, '-', 'Color', 'm', 'LineWidth', 4);
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Geom','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 



load("results_comparison3.mat")

figure(1)

subplot(2,3,3)
hold on
errorbar(mvecs, err_detective./trG, SD_detective./trG, '-', 'Color', 'm', 'LineWidth', 4);
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Gaps','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 
 



load("results_comparison4.mat")

figure(1)

subplot(2,3,4)
hold on
errorbar(mvecs, err_detective./trG, SD_detective./trG, '-', 'Color', 'm', 'LineWidth', 4);
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example RBF','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 



load("results_comparison5.mat")

figure(1)

subplot(2,3,5)
hold on
errorbar(mvecs, err_detective./trG, SD_detective./trG, '-', 'Color', 'm', 'LineWidth', 4);
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(1/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 



load("results_comparison6.mat")

figure(1)

subplot(2,3,6)
hold on
h1 = errorbar(mvecs, err_detective./trG, SD_detective./trG, '-', 'Color', 'm', 'LineWidth', 4);
h2 = errorbar(mvecs+10, err_HS./trG, SD_HS./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
h3 = errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
h4 = errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
h5 = errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '-', 'Color', colors(1,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(3/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


ax_leg = axes('Position',[0.1 0.01 0.8 0.05],'Visible','off'); 
lgd = legend(ax_leg, [h1 h2 h3 h4 h5], {'log-det-ective','half samples','one-sample','LowRank','Plain G-H'}, ...
    'Orientation','horizontal','FontSize',18,'NumColumns',5);
lgd.Box = 'off';  