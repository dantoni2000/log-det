%% This code computes the errors produced by the mixed strategy, for values β = {1/8, 1/4, 1/2, 3/4, 7/8}, with diagonal matrices

clear all
close all
warning off

addpath('tools')

% number of tries
T = 50;

% number of the experiment
c = 1;

% dimension of the matrix
n = 4000; 

% creation of the matrix
[G,Q,trG] = matrices(c,n);
% A = Q*G*Q';
A = G;

% budget of matvecs for the low-rank approximation
l_min = 100; l_max = 1000; step = 100; lgt = 1 + (l_max - l_min) / step;

% iterations of Lanczos
m = 10;

% Preallocation of the matrices for the errors (T x lgt)
ERR_beta18 = zeros(T,lgt);
ERR_beta14 = zeros(T,lgt);
ERR_beta12 = zeros(T,lgt);
ERR_beta34 = zeros(T,lgt);
ERR_beta78 = zeros(T,lgt);

% Final averaged errors
beta18 = zeros(lgt,1); SD_beta18  = zeros(lgt,1);
beta14 = zeros(lgt,1); SD_beta14  = zeros(lgt,1);
beta12 = zeros(lgt,1); SD_beta12  = zeros(lgt,1);
beta34 = zeros(lgt,1); SD_beta34  = zeros(lgt,1);
beta78 = zeros(lgt,1); SD_beta78  = zeros(lgt,1);

mvecs = zeros(lgt,1);

% comparison of the errors changing the total matvecs
for l = l_min:step:l_max

    mvecs(l/step) = l   
    
    for t = 1:T

        % beta = 1/8
        [~,Beta18] = log_det_ective(A,l,m,1/8);
        ERR_beta18(t,l/step) = abs(Beta18 - trG);

        % beta = 1/4
        [~,Beta14] = log_det_ective(A,l,m,1/4);
        ERR_beta14(t,l/step) = abs(Beta14 - trG);

        % beta = 1/2
        [~,Beta12] = log_det_ective(A,l,m,1/2);
        ERR_beta12(t,l/step) = abs(Beta12 - trG);

        % beta = 3/4
        [~,Beta34] = log_det_ective(A,l,m,3/4);
        ERR_beta34(t,l/step) = abs(Beta34 - trG);

        % beta = 7/8
        [~,Beta78] = log_det_ective(A,l,m,7/8);
        ERR_beta78(t,l/step) = abs(Beta78 - trG);

    end

end
 
% computation of the errors (10% trimmed)
trim = floor(0.1*T);
keep = trim+1 : T-trim;

Xs = sort(ERR_beta18,1);
beta18 = mean(Xs(keep,:),1).';
SD_beta18 = std(Xs(keep,:),0,1).';

Xs = sort(ERR_beta14,1);
beta14 = mean(Xs(keep,:),1).';
SD_beta14 = std(Xs(keep,:),0,1).';

Xs = sort(ERR_beta12,1);
beta12 = mean(Xs(keep,:),1).';
SD_beta12 = std(Xs(keep,:),0,1).';

Xs = sort(ERR_beta34,1);
beta34 = mean(Xs(keep,:),1).';
SD_beta34 = std(Xs(keep,:),0,1).';

Xs = sort(ERR_beta78,1);
beta78 = mean(Xs(keep,:),1).';
SD_beta78 = std(Xs(keep,:),0,1).';

subfolder = 'betas_comparison';
filename = fullfile(subfolder, sprintf('results_betas%d.mat', c));
save(filename, ...
    'ERR_beta18','ERR_beta14','ERR_beta12','ERR_beta34','ERR_beta78', ...
    'beta18','SD_beta18', ...
    'beta14','SD_beta14', ...
    'beta12','SD_beta12', ...
    'beta34','SD_beta34', ...
    'beta78','SD_beta78', ...
    'trG');

% % eventual plots
% figure(1)
% semilogy(diag(G), 'LineWidth', 5)
% xlabel('$n$','interpreter','Latex','fontsize',18)
% ylabel('eigenvalues','fontsize',18)
% title('Eigenvalues of the matrix','fontsize',18)
% legend('$\lambda(A)$','interpreter','Latex','fontsize',18)
% ax = gca;
% ax.FontSize = 30; 
% 
% colors = [
%     1.0 0.0 0.0;
%     1.0 0.5 0.0;
%     1.0 1.0 0.0;
%     0.0 0.8 0.0;
%     0.0 0.3 1.0
% ];
% 
% figure(2)
% hold on
% errorbar(mvecs, beta18./trG, SD_beta18./trG, '-', 'Color',_
