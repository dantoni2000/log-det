%% This code computes the stochastic Lanczos quadrature of the matrix log(A+I) via the splitting trace log(A+I) = trace log(P) + trace log(sqrt(P)(A+I)sqrt(P)), using the slq to the residual term

function [its,Ptr] = stochastic_Preconditioned_Lanczos_quadrature(A,l,N,m,tol_beta,tolrel,tolabs)

% input:

% A = input matrix, 
% l = number of matvecs for the Nyström approximation,
% N = number of samples for the Girard-Hutchinson estimator,
% m = number of steps of the Lanczos method,
% tol_beta = tolerance for beta (1e-12 default),
% tolrel = relative tolerance (1e-8 default),
% tolabs = absolute tolerance (1e-8 default).


% output: 

% its = number of total iterations, 
% Ptr = trace log(P) + slq(log(sqrt(P)(A+I)sqrt(P))) ≈ trace log(A+I).


[n,~] = size(A);

Ptr = 0; Pp = 0;
Pit = zeros(N,1);

U = 1; Lhat = 0;

if nargin < 5 || isempty(tol_beta)
    tol_beta = 1e-12;
end
if nargin < 6 || isempty(tolrel)
    tolrel = 1e-8;
end
if nargin < 7 || isempty(tolabs)
    tolabs = 1e-8;
end

% if l>0, compute the Nyström approximation, compute the exact value of trace log(P)
if l>0
    [U,Lhat,~] = Nystrom(A,l);
    Pp = sum(log(diag(Lhat+eye(l,l))),"all");
end
  
%compute the stochastic Lanczos quadrature of the residual matrix for any sample w_i
for i = 1:N
    v = randn(n,1);
    [Pit(i),tr] = Preconditioned_Lanczos_log(A,v,m,U,Lhat,tol_beta,tolrel,tolabs);
    Ptr = Ptr + tr(end-1,1);
end
    
its = 1/N * sum(Pit);
Ptr = 1/N * Ptr + Pp;
