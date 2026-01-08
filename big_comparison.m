%% This code computes the errors produced by the strategies described in the paper, with diagonal matrices

clear all
close all
warning off

addpath('tools')

% number of tries
T = 100;

% number of the experiment 
c = 1;

% iterations of Lanczos
m = 50;

% choice of parameters
if c == 1
    mi = 1e+2;
    d = 5; 
    l = d;
else 
    mi = 1e+2;
    d = 3;
    alpha = 1;
end

n_min = 1000; n_max = 10000; step = 1000; lgt = 1 + (n_max - n_min) / step;
X = randn(d,n_max); 
mvecs = zeros(lgt,1);

% Preallocation of the matrices for the errors of the methods (T x lgt)
ERR_detective = zeros(T,lgt);
ERR_HS        = zeros(T,lgt);
ERR_1S        = zeros(T,lgt);
ERR_LowRank   = zeros(T,lgt);

% Final averaged errors
err_detective = zeros(lgt,1); SD_detective  = zeros(lgt,1);
err_HS        = zeros(lgt,1); SD_HS          = zeros(lgt,1);
err_1S        = zeros(lgt,1); SD_1S          = zeros(lgt,1);
err_LowRank   = zeros(lgt,1); SD_LowRank     = zeros(lgt,1);

for n = n_min:step:n_max
    
    budget = floor(n/10);
    mvecs(n/step) = n + m

    D = pdist2(X(:,1:n)', X(:,1:n)', 'squaredeuclidean'); % n x n

    if c == 1
        A  = exp(-D2 / l^2);
    else 
        A = (1 + alpha*D + (alpha^2/3)*(D.^2)) .* exp(-alpha*D);
        A(1:size(A,1)+1:end) = 1;
    end

    [Q,G,~] = svd(A);
    G = mi*G;
    trG = sum(log(diag(G+eye(n,n))),"all");

    for t = 1:T

        % log-det-ective 
        beta = 3/4; 
        [~,tr_d] = log_det_ective(G,budget,m,beta);
        ERR_detective(t,n/step) = abs(tr_d - trG);

        % half samples
        N = floor(budget/(2*m))+1;
        [~,tr_HS] = stochastic_Preconditioned_Lanczos_quadrature(G,floor(budget/2)+m,N,m);
        ERR_HS(t,n/step) = abs(tr_HS - trG);

        % one sample
        N = 1;
        [~,tr_1S] = stochastic_Preconditioned_Lanczos_quadrature(G,budget,N,m);
        ERR_1S(t,n/step) = abs(tr_1S - trG);

        % low rank
        [~,LhatBig,~] = Nystrom(G,budget+m);
        ERR_LowRank(t,n/step) = abs( ...
            sum(log(diag(LhatBig+eye(budget+m,budget+m))),"all") - trG);


    end
end

% computation of the errors (10% trimmed)
trim = floor(0.1*T);
keep = trim+1 : T-trim;

Xs = sort(ERR_detective,1);
err_detective = mean(Xs,1).';
SD_detective  = std(Xs(keep,:),0,1).';

Xs = sort(ERR_HS,1);
err_HS = mean(Xs,1).';
SD_HS  = std(Xs(keep,:),0,1).';

Xs = sort(ERR_1S,1);
err_1S = mean(Xs,1).';
SD_1S  = std(Xs(keep,:),0,1).';

Xs = sort(ERR_LowRank,1);
err_LowRank = mean(Xs,1).';
SD_LowRank  = std(Xs(keep,:),0,1).';


subfolder = 'big_example';
filename = fullfile(subfolder, sprintf('results_big%d.mat',c));
save(filename, ...
    'ERR_detective','ERR_HS','ERR_1S','ERR_LowRank', ...
    'err_detective','SD_detective', ...
    'err_HS','SD_HS', ...
    'err_1S','SD_1S', ...
    'err_LowRank','SD_LowRank', ...
    'trG');

% eventual plots 

colors = [ 0.36 0.14 0.32; 
           1.0 0.5 0.0; 
           0.0 0.8 0.0;
           1.0 0.0 1.0;
           0.0 0.3 1.0];

figure(2)
hold on
errorbar(mvecs, err_detective./abs(trG), SD_detective./abs(trG), '-', 'Color', 'm', 'LineWidth', 4);
errorbar(mvecs+10, err_HS./abs(trG), SD_HS./abs(trG), '-', 'Color', colors(3,:), 'LineWidth', 4);
errorbar(mvecs-10, err_1S./abs(trG), SD_1S./abs(trG), '-', 'Color', colors(5,:), 'LineWidth', 4);
errorbar(mvecs-10, err_LowRank./abs(trG), SD_LowRank./abs(trG), '-', 'Color', colors(2,:), 'LineWidth', 4);
hold off
set(gca, 'YScale', 'log'); 
xlabel('matvecs','interpreter','Latex','fontsize',18)
ylabel('Relative errors','fontsize',18)
title('Comparison of the errors for the strategies','fontsize',18)
legend('error detective', 'error HalfSamples', 'error 1Sample', ...
       'error LowRank','fontsize',18)
ax = gca;
ax.FontSize = 30;
