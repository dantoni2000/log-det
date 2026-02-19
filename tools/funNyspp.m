%% This code runs the funNys++ algorithm.
%%  The implementation follows from:
% [Algorithm 5.1] D.Persson, D.Kressner - Randomized low-rank approximation of monotone matrix functions 
% https://doi.org/10.1137/22M1523923

% input: 

% A = input matrix, 
% s = number of matvecs for the Nyström approximation,
% N = number of samples for the Girard-Hutchinson estimator,
% m = number of steps of the Lanczos method,
% tolrel = relative tolerance (1e-8 default),
% tolabs = absolute tolerance (1e-8 default).

% output: 

% tr ≈ trace log(A+I)

function [its,Ptr] = funNyspp(A,s,N,m,tolrel,tolabs)

if nargin < 5 || isempty(tolrel)
        tolrel = 1e-8; 
end
if nargin < 6 || isempty(tolabs)
        tolabs = 1e-8; 
end

[n,~] = size(A);

Ptr = 0;
Pit = zeros(N,1);

[U,Lhat] = Nystrom(A,s);
Pp = sum(log(diag(Lhat+eye(s,s))),"all");
    
for i = 1:N
    v = randn(n,1); w = U'*v;
    [Pit(i),tr] = Preconditioned_Lanczos_log(A,v,m,1,0,1e-12,tolrel,tolabs);
    tr_P = w' * (log(diag(Lhat+eye(s,s))).*w);
    Ptr = Ptr + (tr(end-1,1)-tr_P);
end
    
its = 1/N * sum(Pit);
Ptr = 1/N * Ptr + Pp;
