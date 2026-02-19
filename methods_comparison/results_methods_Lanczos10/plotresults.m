mvecs = 110:100:1010;

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
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, 'd-', 'Color', colors(3,:), 'LineWidth', 2);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 2);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '--', 'Color', colors(1,:), 'LineWidth', 2);
errorbar(mvecs, err_kraw_power./trG, SD_kraw_power./trG, '*-', 'Color', 'k', 'LineWidth', 2);
errorbar(mvecs, err_funNyspp./trG, SD_funNyspp./trG, 'o-', 'Color', 'c', 'LineWidth', 2);
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
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, 'd-', 'Color', colors(3,:), 'LineWidth', 2);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 2);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '--', 'Color', colors(1,:), 'LineWidth', 2);


% errorbar(mvecs(1:5), err_kraw_power(1:5)./trG, SD_kraw_power(1:5)./trG, '*-', 'Color', 'k', 'LineWidth', 2);

grayTrans = [0.75 0.75 0.75];   % grigio un po' più scuro

% 1:4 in nero (marker+linea+errorbar)
errorbar(mvecs(1:4), err_kraw_power(1:4)./trG, SD_kraw_power(1:4)./trG, '*-', ...
    'Color', 'k', 'LineWidth', 2);

% linea 4 -> 5 in grigio
plot(mvecs(4:5), err_kraw_power(4:5)./trG, '-', ...
    'Color', grayTrans, 'LineWidth', 2);

% punto 5 in grigio (marker+errorbar, senza linea)
errorbar(mvecs(5), err_kraw_power(5)./trG, SD_kraw_power(5)./trG, '*', ...
    'Color', grayTrans, 'LineWidth', 2, 'LineStyle', 'none');



errorbar(mvecs, err_funNyspp./trG, SD_funNyspp./trG, 'o-', 'Color', 'c', 'LineWidth', 2);
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
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, 'd-', 'Color', colors(3,:), 'LineWidth', 2);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 2);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '--', 'Color', colors(1,:), 'LineWidth', 2);


grayTrans = [0.75 0.75 0.75];   % grigio un po' più scuro

% 1:7 in nero (marker+linea+errorbar)
errorbar(mvecs(1:7), err_kraw_power(1:7)./trG, SD_kraw_power(1:7)./trG, '*-', ...
    'Color', 'k', 'LineWidth', 2);

% linea 7 -> 8 in grigio
plot(mvecs(7:8), err_kraw_power(7:8)./trG, '-', ...
    'Color', grayTrans, 'LineWidth', 2);

% punto 8 in grigio (marker+errorbar, senza linea)
errorbar(mvecs(8), err_kraw_power(8)./trG, SD_kraw_power(8)./trG, '*', ...
    'Color', grayTrans, 'LineWidth', 2, 'LineStyle', 'none');


errorbar(mvecs, err_funNyspp./trG, SD_funNyspp./trG, 'o-', 'Color', 'c', 'LineWidth', 2);
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
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, 'd-', 'Color', colors(3,:), 'LineWidth', 2);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 2);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '--', 'Color', colors(1,:), 'LineWidth', 2);
errorbar(mvecs, err_kraw_power./trG, SD_kraw_power./trG, '*-', 'Color', 'k', 'LineWidth', 2);
errorbar(mvecs, err_funNyspp./trG, SD_funNyspp./trG, 'o-', 'Color', 'c', 'LineWidth', 2);
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
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, 'd-', 'Color', colors(3,:), 'LineWidth', 2);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 2);
errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '--', 'Color', colors(1,:), 'LineWidth', 2);
errorbar(mvecs, err_kraw_power./trG, SD_kraw_power./trG, '*-', 'Color', 'k', 'LineWidth', 2);
errorbar(mvecs, err_funNyspp./trG, SD_funNyspp./trG, 'o-', 'Color', 'c', 'LineWidth', 2);
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
h2 = errorbar(mvecs+10, err_HS./trG, SD_HS./trG, 'd-', 'Color', colors(3,:), 'LineWidth', 2);
h3 = errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
h4 = errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 2);
h5 = errorbar(mvecs(err_PlainGH < 10^8), err_PlainGH(err_PlainGH < 10^8)./trG, SD_PlainGH(err_PlainGH < 10^8)./trG, '--', 'Color', colors(1,:), 'LineWidth', 2);
h6 = errorbar(mvecs, err_kraw_power./trG, SD_kraw_power./trG, '*-', 'Color', 'k', 'LineWidth', 2);
h7 = errorbar(mvecs, err_funNyspp./trG, SD_funNyspp./trG, 'o-', 'Color', 'c', 'LineWidth', 2);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(3/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


ax_leg = axes('Position',[0.1 0.01 0.8 0.05],'Visible','off'); 
lgd = legend(ax_leg, [h1 h2 h3 h4 h5 h6 h7], {'log-det-ective','half samples','one-sample','low-rank','SLQ','Krylov-aware','funNys++'}, ...
    'Orientation','horizontal','FontSize',30,'NumColumns',7);
lgd.Box = 'off';