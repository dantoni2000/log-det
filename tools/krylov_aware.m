%% This code runs the Krylov-aware algorithm.
%%  The implementation follows from:
% [Algorithm 3.1] T.Chen, E.Hallman - Krylov-Aware Stochastic Trace Estimation
% https://doi.org/10.1137/22M1494257

% input: 

% A = input matrix, 
% b = block parameter,
% q = Krylov extra step parameter,
% n = Krylov step parameter,
% m = number of samples for the residual term,
% ort = ortogonalization parameter.

% output: 

% tr ≈ trace log(A+I)

function tr = krylov_aware(A,b,q,n,m,ort)


d = size(A,1);
Omega = randn(d,b); 

if ort == 1 
    Omega = orth(Omega); 
end 

[Q, T] = block_Lanczos(A,Omega,q+1,n);
Q = Q(:,1:(q+1)*b);

[U,S] = eig(T);
logS = U * (diag(log(diag(S)+1)).' * U');
logs = logS(1:(q+1)*b,1:(q+1)*b); 
logs = (logs + logs') / 2;
tr_defl = trace(logs);
% [~,s] = eig(logs);
% tr_defl = sum(diag(s));

Psi = randn(d,m); 

if ort == 1
    Psi = orth(Psi);
end

Y = Psi - Q*(Q'*Psi);

tr_rem = 0;

for i = 1:m
    [~, Ti] = block_Lanczos(A,Y(:,i),1,n);
    [V,L] = eig(Ti);
    li = V*(diag(log(diag(L)+1)).' * V');

    if ort ==1
        c = (d-(q+1)*b);
    else 
        c = norm(Y(:,i))^2;
    end

    tr_rem = tr_rem + (c / m) * li(1,1);
end

tr = tr_defl + tr_rem;

end


%% 

function [Q1,T] = block_Lanczos(A,Z,s,r,varargin)

% Code for the block Lanczos method

if nargin == 5 && varargin{1} == true
    
    orthogonality_loss = zeros(1,s);
    three_term_rr_loss = zeros(1,s);
    
end

q = s + r;

% Find block size
b = size(Z,2);

% First iteration
[Q,R0] = qr(Z,0);
R = R0;

for k = 0:q-1
    
    %Set Qk
    Qk = Q(:,(1+k*b):((k+1)*b));
    
    if k == 0
        
        Z = A*Qk;
        
    else
        
        % Set Qk-1
        Qkm1 = Q(:,(1+(k-1)*b):(k*b));
        
        Z = A*Qk - Qkm1*R';
        
    end
    
    % Obtain diagonal block in T
    M = Qk'*Z;
    
    if k == 0
        
        T = M;
        
    else
        
        m = size(T,1);
        T = [T zeros(m,b);zeros(b,m) M];
        b_old = size(R,2);
        T((end-b+1):end,(end-b-b_old+1):(end-b)) = R;
        T((end-b-b_old+1):(end-b),(end-b+1):end) = R';
        
    end
    
    if k <= s-1
        
        Q1 = Q;
        
    end
    
    if k == q-1
        
        if nargin == 5 && varargin{1} == true
            varargout{1} = orthogonality_loss;
            varargout{2} = three_term_rr_loss;
        end
        return
        
    end
        
    
    % Reorthogonalization
    Z = Z - Qk*M;
    
    % Double reorthgonalization
    if k > 0

        Z = Z - Q(:,1:(k*b))*(Q(:,1:(k*b))'*Z);

    end
    
    % Obtain next block
    [Qkp1,R,P] = qr(Z,0);
    
    % Check if rank deficient
    if min(abs(diag(R))) <= (1e-10)*max(abs(diag(R)))

        warning('Rank deficient block detected')

        % New block size
        b = max(find(abs(diag(R)) > (1e-10)*max(abs(diag(R)))));

        % Truncate
        R = R(:,1:b);
        Qkp1 = Qkp1(:,1:b);

    end
    
    % Permute back
    invP = zeros(1,length(P));
    invP(P) = 1:length(P);
    R = R(:,invP);
    if nargin == 5 && varargin{1} == true
        
        orthogonality_loss(k+1) = norm(Q'*Q-eye(size(Q,2)),'fro');
        Ek = eye(size(Q,2)); 
        Ek = Ek(:,((end-b)+1):end);
        three_term_rr_loss(k+1) = norm(A*Q - Q*T-Qkp1*R*Ek');
        
    end
        
    Q = [Q Qkp1];
    
    
end

end