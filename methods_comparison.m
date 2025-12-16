%% This code computes the errors produced by the strategies described in the paper, with diagonal matrices

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

% creation of the matrix: A = Q*G*Q';
[G, Q, trG] = matrices(c,n);
% A = Q*G*Q';
A = G; 

% budget of matvecs for the low-rank approximation
l_min = 100; l_max = 1000; step = 100; lgt = 1 + (l_max - l_min) / step;

% iterations of Lanczos
m = 5;

% Preallocation of the matrices for the errors of the methods (T x lgt)
ERR_detective = zeros(T,lgt);
ERR_HS        = zeros(T,lgt);
ERR_1S        = zeros(T,lgt);
ERR_LowRank   = zeros(T,lgt);
ERR_funNys    = zeros(T,lgt);
ERR_PlainGH   = zeros(T,lgt);

% Final averaged errors
err_detective = zeros(lgt,1); SD_detective  = zeros(lgt,1);
err_HS        = zeros(lgt,1); SD_HS          = zeros(lgt,1);
err_1S        = zeros(lgt,1); SD_1S          = zeros(lgt,1);
err_LowRank   = zeros(lgt,1); SD_LowRank     = zeros(lgt,1);
err_funNys    = zeros(lgt,1); SD_funNys      = zeros(lgt,1);
err_PlainGH   = zeros(lgt,1); SD_PlainGH     = zeros(lgt,1);

S = [];
mvecs = zeros(lgt,1);

% comparison of the errors changing the total matvecs
for l = l_min:step:l_max

    mvecs(l/step) = l
        
    for t = 1:T

        % log-det-ective 
        beta = 3/4; 
        [~,tr_d] = log_det_ective(A,l,m,beta);
        ERR_detective(t,l/step) = abs(tr_d - trG);

        % half samples
        N = l/(2*m);
        [~,tr_HS] = stochastic_Preconditioned_Lanczos_quadrature(A,l/2+m,N,m);
        ERR_HS(t,l/step) = abs(tr_HS - trG);

        % one sample
        N = 1;
        [~,tr_1S] = stochastic_Preconditioned_Lanczos_quadrature(A,l,N,m);
        ERR_1S(t,l/step) = abs(tr_1S - trG);

        % low rank
        [~,LhatBig,~] = Nystrom_sanity(A,l+m);
        ERR_LowRank(t,l/step) = abs( ...
            sum(log(diag(LhatBig+eye(l+m,l+m))),"all") - trG );

        % funNys++
        N = l/(2*m);
        [~,tr_fN] = funNyspp(A,l/2+m,N,m);
        ERR_funNys(t,l/step) = abs(tr_fN - trG);

        % SLQ
        FrG = sqrt(sum(log(diag(G+eye(n,n))).^2,"all"));
        proxi = m / (l+m) * FrG; 
        m_conv = stochastic_Preconditioned_Lanczos_quadrature(A,0,1,max(l+m,n),[],[],proxi);
        m_conv = max(m,m_conv);
        s = floor((l+m)/m_conv);

        [~,tr_PGH] = stochastic_Preconditioned_Lanczos_quadrature(A,0,s,m_conv);
        ERR_PlainGH(t,l/step) = abs(tr_PGH - trG);

    end

    S = [S;s];
end

% computation of the errors (10% trimmed)
trim = floor(0.1*T);
keep = trim+1 : T-trim;

Xs = sort(ERR_detective,1);
err_detective = mean(Xs(keep,:),1).';
SD_detective  = std(Xs(keep,:),0,1).';

Xs = sort(ERR_HS,1);
err_HS = mean(Xs(keep,:),1).';
SD_HS  = std(Xs(keep,:),0,1).';

Xs = sort(ERR_1S,1);
err_1S = mean(Xs(keep,:),1).';
SD_1S  = std(Xs(keep,:),0,1).';

Xs = sort(ERR_LowRank,1);
err_LowRank = mean(Xs(keep,:),1).';
SD_LowRank  = std(Xs(keep,:),0,1).';

Xs = sort(ERR_funNys,1);
err_funNys = mean(Xs(keep,:),1).';
SD_funNys  = std(Xs(keep,:),0,1).';

Xs = sort(ERR_PlainGH,1);
err_PlainGH = mean(Xs(keep,:),1).';
SD_PlainGH  = std(Xs(keep,:),0,1).';

subfolder = 'methods_comparison';
filename = fullfile(subfolder, sprintf('results_comparison%d.mat',c));
save(filename, ...
    'ERR_detective','ERR_HS','ERR_1S','ERR_LowRank','ERR_funNys','ERR_PlainGH', ...
    'err_detective','SD_detective', ...
    'err_HS','SD_HS', ...
    'err_1S','SD_1S', ...
    'err_LowRank','SD_LowRank', ...
    'err_funNys','SD_funNys', ...
    'err_PlainGH','SD_PlainGH', ...
    'trG');

% eventual plots 
figure(1)
semilogy(diag(G), 'LineWidth', 5)
xlabel('$n$','interpreter','Latex','fontsize',18)
ylabel('Eigenvalues','fontsize',18)
title('Eigenvalues of the input matrix','fontsize',18)
legend('$\lambda(A)$','interpreter','Latex','fontsize',18)
ax = gca;
ax.FontSize = 30; 

figure(2)
hold on
errorbar(mvecs,err_detective./trG,SD_detective./trG,'-g', 'LineWidth', 5) 
errorbar(mvecs,err_HS./trG,SD_HS./trG,'-b', 'LineWidth', 5) 
errorbar(mvecs,err_1S./trG,SD_1S./trG,'-r', 'LineWidth', 5) 
errorbar(mvecs,err_LowRank./trG,SD_LowRank./trG,'-k', 'LineWidth', 5) 
errorbar(mvecs,err_funNys./trG,SD_funNys./trG,'-c', 'LineWidth', 5) 
errorbar(mvecs(S>1),err_PlainGH(S>1)./trG,SD_PlainGH(S>1)./trG, ...
    'Color', [0.9290 0.6940 0.1250], 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'LineWidth', 5)
hold off
set(gca, 'YScale', 'log'); 
xlabel('$\ell$','interpreter','Latex','fontsize',18)
ylabel('Relative errors','fontsize',18)
title('Comparison of the errors for the strategies','fontsize',18)
legend('error detective', 'error HalfSamples', 'error 1Sample', ...
       'error LowRank', 'error funNys++', 'error Plain GH','fontsize',18)
ax = gca;
ax.FontSize = 30;
