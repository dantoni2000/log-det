mvecs = 110:100:1010;

% ===================== COLORS =====================
colors = [
    0.0 0.3 1.0   % blue (Fro)
    0.0 0.8 0.0
    1.0 0.5 0.0   % orange (Tr)
];

% ---- pastel orange -> blue for Mix ----
nMix = 9;
purple = [0.22  0.08  0.32];
blue   = [0.00  0.65  0.80];

mixColors = [ ...
    linspace(purple(1), blue(1), nMix)', ...
    linspace(purple(2), blue(2), nMix)', ...
    linspace(purple(3), blue(3), nMix)' ...
];

pastelFactor = 0.45;
mixColors = mixColors + pastelFactor * (1 - mixColors);

figure(1)

% ===================== 1 =====================
load("optimal_bounds_comparison1.mat")
subplot(2,3,1)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, BestErr_Tr./trG,  '-', 'Color', colors(3,:), 'LineWidth', 6);
for i = 1:9
    semilogy(mvecs, BestErr_Mix(:,i)./trG, ...
        'Color', mixColors(i,:), 'LineWidth', 2);
end
hold off
set(gca,'YScale','log','FontSize',15)
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Alg','fontsize',15)

% ===================== 2 =====================
load("optimal_bounds_comparison2.mat")
subplot(2,3,2)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, BestErr_Tr./trG,  '-', 'Color', colors(3,:), 'LineWidth', 6);
for i = 1:9
    semilogy(mvecs, BestErr_Mix(:,i)./trG, ...
        'Color', mixColors(i,:), 'LineWidth', 2);
end
hold off
set(gca,'YScale','log','FontSize',15)
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Geom','fontsize',15)

% ===================== 3 =====================
load("optimal_bounds_comparison3.mat")
subplot(2,3,3)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, BestErr_Tr./trG,  '-', 'Color', colors(3,:), 'LineWidth', 6);
for i = 1:9
    semilogy(mvecs, BestErr_Mix(:,i)./trG, ...
        'Color', mixColors(i,:), 'LineWidth', 2);
end
hold off
set(gca,'YScale','log','FontSize',15)
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Gaps','fontsize',15)

% ===================== 4 =====================
load("optimal_bounds_comparison4.mat")
subplot(2,3,4)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, BestErr_Tr./trG,  '-', 'Color', colors(3,:), 'LineWidth', 6);
for i = 1:9
    semilogy(mvecs, BestErr_Mix(:,i)./trG, ...
        'Color', mixColors(i,:), 'LineWidth', 2);
end
hold off
set(gca,'YScale','log','FontSize',15)
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example RBF','fontsize',15)

% ===================== 5 =====================
load("optimal_bounds_comparison5.mat")
subplot(2,3,5)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, BestErr_Tr./trG,  '-', 'Color', colors(3,:), 'LineWidth', 6);
for i = 1:9
    semilogy(mvecs, BestErr_Mix(:,i)./trG, ...
        'Color', mixColors(i,:), 'LineWidth', 2);
end
hold off
set(gca,'YScale','log','FontSize',15)
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(1/2)','fontsize',15)

% ===================== 6 =====================
load("optimal_bounds_comparison6.mat")
subplot(2,3,6)
hold on
h1 = semilogy(mvecs, BestErr_Fro./trG, '-', 'Color', colors(1,:), 'LineWidth', 6);
h2 = semilogy(mvecs, BestErr_Tr./trG,  '-', 'Color', colors(3,:), 'LineWidth', 6);
for i = 1:9
    semilogy(mvecs, BestErr_Mix(:,i)./trG, ...
        'Color', mixColors(i,:), 'LineWidth', 2);
end
hold off
set(gca,'YScale','log','FontSize',15)
xlabel('matvecs','fontsize',15)
ylabel('Relative errors','fontsize',15)
title('Example Matérn(3/2)','fontsize',15)

% ===================== LEGEND =====================
ax_leg = axes('Position',[0.1 0.01 0.8 0.05],'Visible','off');
lgd = legend(ax_leg,[h1 h2], ...
    {'idealized one-sample','idealized low-rank'}, ...
    'Orientation','horizontal','FontSize',30,'NumColumns',5);
lgd.Box = 'off';
