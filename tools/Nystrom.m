%% This code produces a Nyström approximation of a diagonal matrix A, computes the leave-one-out estimate and check if the approximation is sufficiently good for the one-sample strategy
%% The implementation follows from:
 
% [Algorithm 4.1] E.Epperly, J.Tropp - Efficient error and variance estimation for randomized matrix computations 
% https://doi.org/10.1137/23M1558537

function [V, Lambda, err] = Nystrom(A,l,beta,m)

% input: 
 
% A = diagonal input matrix,
% l = total matvecs with the matrix A, 
% beta = parameter of the mixed strategy (1 default), 
% m = number of steps of Lanczos for the strategies.


% output: 
 
% Nys(A)_l = V Lambda V' spectral decomposition of Nys(A)_l, 
% err ≈ ||A-Nys(A)_l||_F (!) ONLY WHEN beta < 1 (!) 


% Nyström without sanity check: beta = 1, m = 1

if nargin < 3 || isempty(beta)
        beta = 1; 
end
if nargin < 4 || isempty(m)
        m = 1; 
end

d = diag(A);
lb = floor(beta*l);
err = 0;

% Computes the Nyström approximation Nys(A)_βl = V Lambda V'
Omega = randn(length(d), lb);
Y = d .* Omega; 
nu = eps() * norm(Y);
Y = Y + nu * Omega;
[Q, R] = qr(Y, 0);  
H = Omega' * Y;
C = chol((H + H') / 2);
[U, Sigma, ~] = svd(R / C, 'econ');  
Lambda = max(Sigma.^2 - nu * eye(lb), 0);
V = Q * U;
diagHinv = diag(inv(H))';

if beta < 1
    lbb = floor(beta*lb);
    diagHinv_half = diag(inv(H(1:floor(lbb),1:floor(lbb))))';

    % leave-one-out estimate of Nys(A)_βl and Nys(A)_βl/2
    err = norm(((R / C) / (C')) ./ diagHinv, 'fro') / sqrt(lb);
    err_half = norm(((R(1:floor(lbb),1:floor(lbb)) / C(1:floor(lbb),1:floor(lbb))) / (C(1:floor(lbb),1:floor(lbb))')) ./ diagHinv_half, 'fro') / sqrt(lbb);

    % sanity check of the approximation: if Nystrom is working fine, conclude the approximation
    if sqrt(m / ((1-beta)*lb+m)) * err_half >= err
        % m / ((1-beta)*lb+m) * log(1+err_half^2) >= log(1+err^2)
        OOmega = randn(length(d), l-lb);
        YY = d .* OOmega;
        YY = YY + nu * OOmega;
    
        [QQ, RR] = qr([Y YY], 0);
        HH = [H Omega'*YY;
              OOmega'* Y OOmega' * YY];
        CC = chol((HH + HH') / 2);
    
        [U, Sigma, ~] = svd(RR / CC, 'econ');
        
        Lambda = max(Sigma.^2 - nu * eye(size(Sigma)), 0);
        V = QQ * U;
    end

end