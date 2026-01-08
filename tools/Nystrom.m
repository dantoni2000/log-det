%% This code produces a Nyström approximation of a diagonal matrix A, computes the leave-one-out estimate and check if the approximation is sufficiently good for the one-sample strategy
%% The implementation follows from: 
% [Algorithm 4.1] E.Epperly, J.Tropp - Efficient error and variance estimation for randomized matrix computations 
% https://doi.org/10.1137/23M1558537

function [V, Lambda, err] = Nystrom(A, l, beta, m)

% input: 

% A = matrix (dense, sparse, spdiags) or function handle Afun(v)=A*v
% l = total matvecs with the matrix A
% beta = parameter of the mixed strategy (1 default)
% m = number of steps of Lanczos for the strategies.

% output: 

% Nys(A)_l = V Lambda V' spectral decomposition of Nys(A)_l
% err ≈ ||A-Nys(A)_l||_F (!) ONLY WHEN beta < 1 (!)

if nargin < 3 || isempty(beta)
    beta = 1; 
end
if nargin < 4 || isempty(m)
    m = 1; 
end

% Wrap AX as function handle
Afun = @(X) A * X;
n = size(A,1);

lb = floor(beta*l);
err = 0;

% Computes the Nyström approximation Nys(A)_βl = V Lambda V'
Omega = randn(n, lb);
Y = Afun(Omega);

nu = eps * norm(Y,'fro');
Y = Y + nu * Omega;

[Q, R] = qr(Y, 0);
H = Omega' * Y;

C = chol((H + H')/2, 'lower');
[U, Sigma, ~] = svd(R / C', 'econ');

Lambda = max(Sigma.^2 - nu * eye(size(Sigma)), 0);
V = Q * U;

diagHinv = diag(inv(H))';

% Leave-one-out estimator
if beta < 1
    lbb = floor(beta * lb);
    Hh = H(1:lbb,1:lbb);
    diagHinv_half = diag(inv(Hh))';

    err = norm(((R / C') ./ diagHinv), 'fro') / sqrt(lb);
    err_half = norm( ((R(1:lbb,1:lbb) / C(1:lbb,1:lbb)') ./ diagHinv_half), 'fro') / sqrt(lbb);

    % sanity check of the approximation: if Nystrom is working fine, conclude the approximation
    if sqrt(m / ((1-beta)*lb + m)) * err_half >= err
        OOmega = randn(n, l-lb);
        YY = Afun(OOmega);
        YY = YY + nu * OOmega;

        [QQ, RR] = qr([Y YY], 0);

        HH = [H          Omega' * YY;
              OOmega' * Y  OOmega' * YY];

        CC = chol((HH + HH')/2, 'lower');
        [U, Sigma, ~] = svd(RR / CC', 'econ');

        Lambda = max(Sigma.^2 - nu * eye(size(Sigma)), 0);
        V = QQ * U;
    end
end
