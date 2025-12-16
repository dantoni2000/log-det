%% This code computes the bounds for the one-sample, for β-strategies and for the low-rank strategy, assuming optimal preconditioner

clear all
close all
warning off

addpath('tools')

% number of the experiment
c = 6;

% dimension of the matrix
n = 4000;

% creation of the matrix
[G, ~, trG] = matrices(c,n);
% A = Q*G*Q';

% budget of matvecs for the low-rank approximation
l_min = 100; l_max = 1000; step = 100; lgt = 1 + (l_max - l_min) / step;

% iterations of Lanczos
m = 5;
 
% computation of the bounds
for l =l_min:step:l_max

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
