mvecs = 100:100:1000;

load("results_Lanczos1.mat")
colors = [
1.0 0.5 0.0; 
0.0 0.3 1.0;
0.36 0.14 0.32;
];

figure(1)

subplot(2,3,1)
hold on
semilogy(mvecs, Lanczos_10./trG, '-o', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, Lanczos_50./trG, '-o', 'Color', colors(2,:), 'LineWidth', 4);
semilogy(mvecs, traceest./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('$\ell$','Interpreter','Latex','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Alg','fontsize',15)
% legend boxoff
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 


load("results_Lanczos2.mat")

figure(1)

subplot(2,3,2)
hold on
semilogy(mvecs, Lanczos_10./trG, '-o', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, Lanczos_50./trG, '-o', 'Color', colors(2,:), 'LineWidth', 4);
semilogy(mvecs, traceest./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('$\ell$','Interpreter','Latex','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Geom','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 




load("results_Lanczos3.mat")


figure(1)

subplot(2,3,3)
hold on
semilogy(mvecs, Lanczos_10./trG, '-o', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, Lanczos_50./trG, '-o', 'Color', colors(2,:), 'LineWidth', 4);
semilogy(mvecs, traceest./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('$\ell$','Interpreter','Latex','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Gaps','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 

 



load("results_Lanczos4.mat")

figure(1)

subplot(2,3,4)
hold on
semilogy(mvecs, Lanczos_10./trG, '-o', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, Lanczos_50./trG, '-o', 'Color', colors(2,:), 'LineWidth', 4);
semilogy(mvecs, traceest./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('$\ell$','Interpreter','Latex','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example RBF','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 
 



load("results_Lanczos5.mat")

figure(1)

subplot(2,3,5)
hold on
semilogy(mvecs, Lanczos_10./trG, '-o', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, Lanczos_50./trG, '-o', 'Color', colors(2,:), 'LineWidth', 4);
semilogy(mvecs, traceest./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('$\ell$','Interpreter','Latex','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(1/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15; 




load("results_Lanczos6.mat")


figure(1)

subplot(2,3,6)
hold on
h1 = semilogy(mvecs, Lanczos_10./trG, '-o', 'Color', colors(1,:), 'LineWidth', 6);
h2 = semilogy(mvecs, Lanczos_50./trG, '-o', 'Color', colors(2,:), 'LineWidth', 4);
h3 = semilogy(mvecs, traceest./trG, '-o', 'Color', colors(3,:), 'LineWidth', 4);
hold off

set(gca, 'YScale', 'log'); 
xlabel('$\ell$','Interpreter','Latex','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(3/2)','fontsize',15)
ax = gca;             % 'gca' prende l'asse corrente
ax.FontSize = 15;

ax_leg = axes('Position',[0.1 0.01 0.8 0.05],'Visible','off'); 
lgd = legend(ax_leg, [h1 h2 h3], {'Lanc($\hat M$,10)','Lanc($\hat M$,50)','error of the trace estimator'}, ...
    'Orientation','horizontal','interpreter','latex' ,'fontsize',30,'NumColumns',5);
lgd.Box = 'off';