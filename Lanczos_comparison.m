%% This code computes the errors produced by the strategies described in the paper, with diagonal matrices

clear all
close all
% warning off

addpath('tools')

% number of tries
T = 50;

% number of the experiment
c = 1;

% dimension of the matrix
n = 4000; 

% creation of the matrix: A = Q*G*Q';
[G, Q, trG] = matrices(c,n);
% A = Q*G*Q';
A = G; 

% budget of matvecs for the low-rank approximation
l_min = 100; l_max = 1000; step = 100; lgt = 1 + (l_max - l_min) / step;

% iterations of Lanczos
m1 = 10; 
m2 = 50;

% Preallocation of the matrices for the errors (T x lgt)
ERR_Lanczos_10 = zeros(T,lgt);
ERR_Lanczos_50 = zeros(T,lgt);
ERR_traceest    = zeros(T,lgt);

% Final averaged errors
Lanczos_10 = zeros(lgt,1); 
Lanczos_50 = zeros(lgt,1); 
traceest   = zeros(lgt,1);

mvecs = zeros(lgt,1);

% comparison of the errors changing the total matvecs
for l = l_min:step:l_max

    mvecs(l/step) = l;
        
    [V, L, ~] = Nystrom_sanity(A,l);
    invsqrtL = diag(diag(L + eye(size(L))).^-0.5);
    
    % sqrt(P)^-1 = V * sqrt(L+I)^-1 V' + (I - V * V')
    invsqrtP = V * invsqrtL * V' + eye(n) - V*V';
    
    % sqrt(P)^-1 (A+I) sqrt(P)^-1 = sqrt(P)^-1 (A-A_l) sqrt(P)^-1 + I
    PAP = invsqrtP * (A - V*L*V') * invsqrtP;
    [~,Sigma] = svd(PAP);
    tracelogPAP = sum(log(diag(Sigma+eye(n,n))),"all");

    for t = 1:T
        v = randn(n,1);

        true_value = v'*(log(diag(Sigma)+1).*v); 

        ERR_traceest(t,l/step) = abs(tracelogPAP - true_value);

        [~,tr10] = Preconditioned_Lanczos_log(Sigma,v,m1);
        ERR_Lanczos_10(t,l/step) = abs(true_value - tr10(end,1));

        [~,tr50] = Preconditioned_Lanczos_log(Sigma,v,m2);
        ERR_Lanczos_50(t,l/step) = abs(true_value - tr50(end,1));
    end

end

% computation of the errors 

Lanczos_10 = mean(ERR_Lanczos_10,1).';

Lanczos_50 = mean(ERR_Lanczos_50,1).';

traceest = mean(ERR_traceest,1).';

subfolder = 'Lanczos_comparison';
filename = fullfile(subfolder, sprintf('results_Lanczos%d.mat',c));
save(filename, ...
    'ERR_Lanczos_10','ERR_Lanczos_50','ERR_traceest', ...
    'Lanczos_10','Lanczos_50','traceest','trG');

figure(1)
semilogy(mvecs, Lanczos_10)
hold on
semilogy(mvecs, Lanczos_50)
hold on
semilogy(mvecs, traceest)
