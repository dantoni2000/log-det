%% This code runs the Lanczos algorithm for the function log(1+x). This version works for a diagonal matrix A, and general P = U L U^T + I  

function [its,tr]=Preconditioned_Lanczos_log(A,x,m,U,L,tol_beta,tolrel,tolabs)

% input:

% A = diagonal input matrix,
% x = starting vector for the iteration,
% m = number of steps of Lanczos,
% [U, L] such that P = U L U^T + I is the preconditioner
% tol_beta = tolerance for beta (1e-12 default),
% tolrel = relative tolerance (1e-8 default),
% tolabs = absolute tolerance (1e-8 default).


% output: 

% its = number of total iterations, 
% tr ≈ trace log(sqrt(P)(A+I)sqrt(P)).


% Preconditioner defaults
if nargin < 4 || isempty(U)
    U = 1;
end
if nargin < 5 || isempty(L)
    L = 0;
end

% Tolerance defaults
if nargin < 6 || isempty(tol_beta)
    tol_beta = 1e-12;
end
if nargin < 7 || isempty(tolrel)
    tolrel = 1e-8;
end
if nargin < 8 || isempty(tolabs)
    tolabs = 1e-8;
end


its = m;
tr(1,1) = 0; 

% vectorize A
d = diag(A);

% sqrt(L+I)
L = diag(diag(L + eye(size(L))).^-0.5);

% first iteration, exploiting the diagonal structure of A
u(:,1) = x/norm(x);
Pu(:,1) = U*(diag(L).*(U'*u(:,1))) + (u(:,1) - U*(U'*u(:,1)));

APu(:,1) = d .* Pu(:,1) + Pu(:,1); % A is diagonal (!)

w(:,1) = U*(diag(L).*(U'*APu(:,1))) + (APu(:,1) - U*(U'*APu(:,1)));
alpha(1) = u(:,1)'*w(:,1);
r(:,1) = w(:,1) - alpha(1)*u(:,1);
beta(1) = norm(r(:,1));

if beta(1) < tol_beta
    tr(1,1) = norm(x)^2 * log(alpha(1));
    its = 1;
    return
end

u(:,2) = r(:,1)/beta(1);
T(1,1) = alpha(1);

for i=2:m
    oldtr = tr(i-1,1);
    
    Pu(:,i) = U*(diag(L).*(U'*u(:,i))) + (u(:,i) - U*(U'*u(:,i)));
    APu(:,i) = d .* Pu(:,i) + Pu(:,i); % A is diagonal (!)
    w(:,i) = U*(diag(L).*(U'*APu(:,i))) + (APu(:,i) - U*(U'*APu(:,i)));
    
    alpha(i,1) = u(:,i)'*w(:,i);
    r(:,i) = w(:,i) - alpha(i)*u(:,i) - beta(i-1)*u(:,i-1);

    % G-S orthogonalization
    for j=1:i
        r(:,i) = r(:,i) - (u(:,j)'*r(:,i))*u(:,j);
    end

    % G-S reorthogonalization for numerical stability
    for j=1:i
        r(:,i) = r(:,i) - (u(:,j)'*r(:,i))*u(:,j);
    end

    beta(i,1) = norm(r(:,i));
    
    if beta(i) < tol_beta
        T(i,i) = alpha(i); 
        T(i,i-1) = beta(i-1); 
        T(i-1,i) = beta(i-1);
        fT = logm(T(1:i, 1:i));
        tr(i,1) = norm(x)^2 * fT(1,1);
        its = i;
        break
    end
    
    u(:,i+1) = r(:,i)/beta(i);
    T(i,i) = alpha(i,1); 
    T(i,i-1) = beta(i-1,1); 
    T(i-1,i) = beta(i-1,1);
    fT = logm(T);
    tr(i,1) = norm(x)^2 * fT(1,1);

    if abs(tr(i,1) - oldtr) <= tolrel*abs(tr(i,1))
        its = i-1;
        break
    end

    if abs(tr(i,1) - oldtr) <= tolabs
        its = i-1;
        break
    end
end
