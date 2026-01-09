mvecs = 1000:1000:10000

figure(1)

% ===================== 1 =====================
load("results_big1.mat")
colors = [ 0.36 0.14 0.32; 
           1.0 0.5 0.0; 
           0.0 0.8 0.0;
           1.0 0.0 1.0;
           0.0 0.3 1.0];

figure(1)
set(gcf,'Units','centimeters','Position',[2 2 24 11])
subplot(1,2,1)
pos = get(gca,'Position');
pos(2) = pos(2) + 0.1;
pos(4) = pos(4) - 0.1;
set(gca,'Position',pos)
hold on
errorbar(mvecs, err_detective./trG, SD_detective./trG, '-', 'Color', 'm', 'LineWidth', 4);
errorbar(mvecs+10, err_HS./trG, SD_HS./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('$n$','fontsize',15, 'Interpreter','latex')
ylabel('Relative errors','fontsize',15)
title('Example RBF','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


load("results_big2.mat")

figure(1)

subplot(1,2,2)
pos = get(gca,'Position');
pos(2) = pos(2) + 0.1;
pos(4) = pos(4) - 0.1;
set(gca,'Position',pos)
hold on
h1 = errorbar(mvecs, err_detective./trG, SD_detective./trG, '-', 'Color', 'm', 'LineWidth', 4);
h2 = errorbar(mvecs+10, err_HS./trG, SD_HS./trG, '-', 'Color', colors(3,:), 'LineWidth', 4);
h3 = errorbar(mvecs-10, err_1S./trG, SD_1S./trG, '-', 'Color', colors(5,:), 'LineWidth', 4);
h4 = errorbar(mvecs, err_LowRank./trG, SD_LowRank./trG, '-', 'Color', colors(2,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('$n$','fontsize',15, 'Interpreter','latex')
ylabel('Relative errors','fontsize',15)
title('Example Matérn(5/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


ax_leg = axes('Position',[0.1 0.01 0.8 0.05],'Visible','off'); 
lgd = legend([h1 h2 h3 h4], ...
    {'log-det-ective','half samples','one-sample','low-rank'}, ...
    'Orientation','horizontal', 'FontSize',30);  % tieni font grande
lgd.Units = 'normalized';
lgd.Position = [0.12 0.02 0.76 0.06];  % [x y w h], y=0.02 => margine sotto
lgd.Box = 'off';