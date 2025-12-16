%% This code computes the stochastic Lanczos quadrature of the (diagonal) matrix log(A+I) using the lod-det-ective algorithm

function [its,Ptr] = log_det_ective(A,l,m,beta)

% input: 

% A = diagonal input matrix, 
% l = number of matvecs for the Nyström approximation, 
% m = number of steps of the Lanczos method, 
% beta = parameter of the mixed strategy.


% output: 
% its = number of total iterations, 
% Ptr = trace log(P) + slq_beta(log(sqrt(P)(A+I)sqrt(P))) ≈ trace log(A+I).


[n,~] = size(A);
Ptr = 0;

% Compute the Nystrom approximation with βl matvecs
[U,Lhat] = Nystrom_sanity(A,l,beta,m);
[z,~] = size(Lhat);
Pp = sum(log(diag(Lhat+eye(z,z))),"all");

if z>floor(beta*l) %all the budget has been allocated for Nystrom: 1 Sample strategy

    v = randn(n,1);
    [Pit,tr] = Preconditioned_Lanczos_log(A,v,m,U,Lhat);
    Ptr = Ptr + tr(end,1);
    its = sum(Pit);
    Ptr = Ptr + Pp;
    
else % mixed strategy

    N = 1 + floor(l*(1-beta)/m);
    Pit = zeros(N,1);
    for i = 1:N
        v = randn(n,1);
        [Pit(i),tr] = Preconditioned_Lanczos_log(A,v,m,U,Lhat);
        Ptr = Ptr + tr(end,1);
    end
        its = 1/(N)*(sum(Pit));
        Ptr = 1/(N)*(Ptr) + Pp;
        Ptr = Ptr(end);
         
end


