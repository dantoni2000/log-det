%% This code computes the bounds for the one-sample and for the low-rank strategy, assuming exact Lanczos

clear all
close all
warning off

% number of the experiment
c = 1;

% dimension of the matrix
n = 4000;

% creation of the matrix
[G, Q, trG] = matrices(c,n);
% A = Q*G*Q';

% budget of matvecs for the low-rank approximation
l_min = 100; l_max = 1000; step = 100; lgt = 1 + (l_max - l_min) / step;

% iterations of Lanczos
m = 50;


% computation of the bounds: takes the minimum over all possible choices of p and k
for l =l_min:step:l_max
    l + m
    
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
