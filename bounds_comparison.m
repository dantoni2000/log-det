%% This code computes the bounds for the one-sample and for the low-rank strategy, assuming exact Lanczos

clear all
close all
warning off

addpath('tools')

% number of the experiment
c = 3;

% dimension of the matrix
n = 4000;

% creation of the matrix
[G, Q, trG] = matrices(c,n);
% A = Q*G*Q';

% budget of matvecs for the low-rank approximation
l_min = 100; l_max = 1000; step = 100; lgt = 1 + (l_max - l_min) / step;

% iterations of Lanczos
m = 10;


% computation of the bounds: takes the minimum over all possible choices of p and k
for l =l_min:step:l_max
    mvecs(l/step) = l + m
    
    for s = 2:2:l-2
        k = s;
        p = l - k;
     
        squareboundFro(:,s/2) = sqrt(2*(k+p-1)/(p-1)) * sum(log(1+diag(G(k+1:n,k+1:n))))^0.5 *( log(1+ (1+2*k/(p-1))* G(k+1,k+1) + 2*exp(2)*(k+p)/(p^2-1) * sum(diag(G(k+1:n,k+1:n)))) )^0.5 ;
     
    end


    for s = 2:2:l+m-2
        k = s;
        p = l + m - k;
       
        squareboundTr(:,s/2) = (k+p-1)/(p-1) * sum(log(1+diag(G(k+1:n,k+1:n))));

     end

    BestErr_Tr(l/step) = min(squareboundTr');
    BestErr_Fro(l/step) = min(squareboundFro');


end


subfolder = 'bounds_comparison';
filename = fullfile(subfolder, sprintf('bounds_comparison%d.mat',c));
save(filename, 'BestErr_Tr' ,'BestErr_Fro', 'trG');


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
1.0 0.5 0.0; 
0.0 0.8 0.0;
0.0 0.3 1.0
];


figure(2)
hold on
semilogy(mvecs, BestErr_Fro./trG, '-o', 'Color', colors(1,:), 'LineWidth', 5);
semilogy(mvecs, BestErr_Tr./trG, '-o', 'Color', colors(3,:), 'LineWidth', 5);
hold off

set(gca, 'YScale', 'log'); 
xlabel('matvecs','fontsize',18)
ylabel('Relative errors','fontsize',18)
title('Example Alg','fontsize',18)
legend('one-sample', 'low-rank', 'fontsize',18)
ax = gca;
ax.FontSize = 30;