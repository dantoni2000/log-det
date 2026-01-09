%% This code computes the bounds for the one-sample, for β-strategies and for the low-rank strategy, assuming optimal preconditioner

clear all
close all
warning off

addpath('tools')

% number of the experiment
c = 3;

% dimension of the matrix
n = 4000;

% creation of the matrix
[G, ~, trG] = matrices(c,n);
% A = Q*G*Q';

% budget of matvecs for the low-rank approximation
l_min = 100; l_max = 1000; step = 100; lgt = 1 + (l_max - l_min) / step;

% iterations of Lanczos
m = 10;
 
% computation of the bounds
for l =l_min:step:l_max
        mvecs(l/step) = l+m
        
        BestErr_Tr(l/step) = sum(log(1+diag(G(l+m+1:n,l+m+1:n))));

        BestErr_Fro(l/step) = sqrt(2* sum( log(1+diag(G(l+1:n,l+1:n))).^2 ) );

        for alpha = 1:9
            
            beta = 0.1*alpha;
            lb = l*beta;
            BestErr_Mix(l/step,alpha) = sqrt(2*m/(l-lb+m) * sum( log(1+diag(G(lb+1:n,lb+1:n))).^2 ) );

        end

end

subfolder = 'bounds_comparison';
filename = fullfile(subfolder, sprintf('optimal_bounds_comparison%d.mat',c));
save(filename, 'BestErr_Tr' ,'BestErr_Fro', 'BestErr_Mix', 'trG');


% eventual plots 
figure(1)
semilogy(diag(G), 'LineWidth', 5)
xlabel('$n$','interpreter','Latex','fontsize',18)
ylabel('Eigenvalues','fontsize',18)
title('Eigenvalues of the input matrix','fontsize',18)
legend('$\lambda(A)$','interpreter','Latex','fontsize',18)
ax = gca;
ax.FontSize = 30;


colors = [
    0.0 0.3 1.0   % blue (Fro)
    0.0 0.8 0.0
    1.0 0.5 0.0   % orange (Tr)
];

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

figure(2)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-', 'Color', colors(1,:), 'LineWidth', 6);
semilogy(mvecs, BestErr_Tr./trG,  '-', 'Color', colors(3,:), 'LineWidth', 6);
for i = 1:9
    semilogy(mvecs, BestErr_Mix(:,i)./trG, ...
        'Color', mixColors(i,:), 'LineWidth', 2);
end
hold off
set(gca,'YScale','log','FontSize',18)
xlabel('matvecs','fontsize',18)
ylabel('Relative errors','fontsize',18)
title('Example Alg','fontsize',18)
ax = gca;
ax.FontSize = 30;