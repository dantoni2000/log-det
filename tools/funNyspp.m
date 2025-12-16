%% This code computes the stochastic Lanczos quadrature of the (diagonal) matrix log(A+I) using the funNys++ algorithm

function [its,Ptr] = funNyspp(A,l,N,m)

% input: 
 
% A = diagonal input matrix, 
% l = number of matvecs for the Nyström approximation, 
% N = number of samples for the Girard-Hutchinson estimator, 
% m = number of steps of the Lanczos method.


% output: 

% its = number of total iterations, 
% Ptr = trace log (P) + slq( log(A+I) - log(P) ) ≈ trace log(A+I).


[n,~] = size(A);

Ptr = 0;
Pit = zeros(N,1);

[U,Lhat,~] = Nystrom_sanity(A,l);
Pp = sum(log(diag(Lhat+eye(l,l))),"all");
    
for i = 1:N
    v = randn(n,1); w = U'*v;
    [Pit(i),tr] = Preconditioned_Lanczos_log(A,v,m);
    tr_P = w' * (log(diag(Lhat+eye(l,l))).*w);
    Ptr = Ptr + (tr(end,1)-tr_P);
end
    
its = 1/N * sum(Pit);
Ptr = 1/N * Ptr + Pp;
