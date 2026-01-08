%% This code runs the Preconditioned Lanczos approximation of trace log( sqrt(P) * (A+I) * sqrt(P) ). % Efficient, generic version: works for dense, sparse, diagonal, or spdiags matrices

function [its, tr] = Preconditioned_Lanczos_log(A, x, m, U, L, tol_beta, tolrel, tolabs)

% Inputs:
%   A        = matrix (dense, sparse, spdiags) or function handle Afun(v)=A*v
%   x        = starting vector
%   m        = number of Lanczos steps
%   U, L     = preconditioner P = U*L*U' + I, L diagonal
%   tol_beta = tolerance for beta (default 1e-12)
%   tolrel   = relative tolerance (default 1e-8)
%   tolabs   = absolute tolerance (default 1e-8)

% Outputs:
%   its = number of Lanczos iterations performed
%   tr ≈ trace log(sqrt(P)(A+I)sqrt(P)).

% Defaults
if nargin < 4 || isempty(U), U = 1; end
if nargin < 5 || isempty(L), L = 0; end
if nargin < 6 || isempty(tol_beta), tol_beta = 1e-12; end
if nargin < 7 || isempty(tolrel), tolrel = 1e-6; end
if nargin < 8 || isempty(tolabs), tolabs = 1e-6; end

its = m;
tr(1,1) = 0;

% Wrap Ax as function handle
if isdiag(A)
    d = diag(A);
    Afun = @(v) d .* v + v;
else
    Afun = @(v) A*v + v;
end

% Preconditioner function, optimized
L_diag = 1 ./ sqrt(diag(L) + 1); % sqrt(L+I)^-0.5
Pu_fun = @(v) local_Pu(v, U, L_diag);

% Initialize Lanczos
u(:,1) = x / norm(x);
Pu = Pu_fun(u(:,1));
APu = Afun(Pu);
w = Pu_fun(APu);

alpha(1) = u(:,1)' * w;
r = w - alpha(1)*u(:,1);
beta(1) = norm(r);

if beta(1) < tol_beta
    tr(1,1) = norm(x)^2 * log(alpha(1));
    its = 1;
    return
end

u(:,2) = r / beta(1);
T = zeros(m,m);
T(1,1) = alpha(1);

% Lanczos loop
for i = 2:m
    oldtr = tr(i-1,1);

    Pu = Pu_fun(u(:,i));
    APu = Afun(Pu);
    w = Pu_fun(APu);

    alpha(i) = u(:,i)'*w;
    r = w - alpha(i)*u(:,i) - beta(i-1)*u(:,i-1);

    % Gram-Schmidt + reorthogonalization
    for j = 1:i
        r = r - (u(:,j)'*r)*u(:,j);
    end
    for j = 1:i
        r = r - (u(:,j)'*r)*u(:,j);
    end

    beta(i) = norm(r);

    % Breakdown
    if beta(i) < tol_beta
        T(i,i) = alpha(i);
        T(i,i-1) = beta(i-1);
        T(i-1,i) = beta(i-1);
        fT = logm(T(1:i,1:i));
        tr(i,1) = norm(x)^2 * fT(1,1);
        its = i;
        break
    end

    u(:,i+1) = r / beta(i);

    % Update T and tr
    T(i,i) = alpha(i);
    T(i,i-1) = beta(i-1);
    T(i-1,i) = beta(i-1);
    fT = logm(T(1:i,1:i));
    tr(i,1) = norm(x)^2 * fT(1,1);

    % Convergence check
    if abs(tr(i,1) - oldtr) <= tolrel*abs(tr(i,1)) || abs(tr(i,1) - oldtr) <= tolabs
        its = i-1;
        break
    end
end

end

%% Helper function: optimized Pu
function y = local_Pu(v, U, L_diag)
tmp = U' * v;                 % compute once
y = U * (L_diag .* tmp) + (v - U * tmp);
end