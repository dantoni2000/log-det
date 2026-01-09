mvecs = 110:100:1010;

load("bounds_comparison1.mat")
colors = [
1.0 0.5 0.0; 
0.0 0.8 0.0;
0.0 0.3 1.0
];

figure(1)

subplot(2,3,1)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-o', 'Color', colors(1,:), 'LineWidth', 4);
semilogy(mvecs, BestErr_Tr./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Alg','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


load("bounds_comparison2.mat")

figure(1)

subplot(2,3,2)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-o', 'Color', colors(1,:), 'LineWidth', 4);
semilogy(mvecs, BestErr_Tr./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Geom','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 




load("bounds_comparison3.mat")


figure(1)

subplot(2,3,3)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-o', 'Color', colors(1,:), 'LineWidth', 4);
semilogy(mvecs, BestErr_Tr./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Gaps','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 

 



load("bounds_comparison4.mat")

figure(1)

subplot(2,3,4)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-o', 'Color', colors(1,:), 'LineWidth', 4);
semilogy(mvecs, BestErr_Tr./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example RBF','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 
 



load("bounds_comparison5.mat")

figure(1)

subplot(2,3,5)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-o', 'Color', colors(1,:), 'LineWidth', 4);
semilogy(mvecs, BestErr_Tr./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(1/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 




load("bounds_comparison6.mat")


figure(1)

subplot(2,3,6)
hold on
h1 = semilogy(mvecs, BestErr_Fro./trG, '-o', 'Color', colors(1,:), 'LineWidth', 4);
h2 = semilogy(mvecs, BestErr_Tr./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(3/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 

ax_leg = axes('Position',[0.1 0.01 0.8 0.05],'Visible','off'); 
lgd = legend(ax_leg, [h1 h2], {'one-sample','LowRank'}, ...
    'Orientation','horizontal','FontSize',18,'NumColumns',5);
lgd.Box = 'off';