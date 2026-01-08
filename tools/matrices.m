%% This code produces the spectral decomposition and the exact tracelog of the matrices described in Section 5. 

function [G, Q, trG] = matrices(c,n)

% input: 

% c = number of the example, 
% n = dimension of the matrix.


% output: 
 
% A = Q G Q' spectral decomposition of A, 
% trG = trace log(A+I).


if c == 1 %Example Alg

    Q = randn(n,n); [Q,~] = qr(Q);
    g = linspace(1,n,n)';
    mi = 10^2;
    G = mi*diag(1./(g).^2);


elseif c == 2 %Example Geom

    Q = randn(n,n); [Q,~] = qr(Q);
    g = linspace(1,n,n)';
    mi = 10^6;
    G = mi*diag(exp(-.1*g));


elseif c == 3 %Example Gaps

   A = zeros(n,n); 
   l = 1e-2; h = 1e+2; k = 1e+0; p = 1e-6; 
   mi = 10^6;

    for j=1:200
        x = sprand(n,1,0.01);
        A = A + (h/j^2*x) * x';
    end

    for j=201:400
        x = sprand(n,1,0.01);
        A = A + (k/j^2*x) * x';
    end

    for j=401:600
        x = sprand(n,1,0.01);
        A = A + (l/j^2*x) * x';
    end

    for j=601:4000
        x = sprand(n,1,0.01);
        A = A + (p/j^2*x) * x';
    end

    [Q,G,~]=svd(A);
    G = mi*G;

        
elseif c == 4 %Example RBF

    x = randn(1,n);
    l = .01;
    mi = 10^2;
    for i = 1:n
        for j = 1:n
            A(i,j) = exp(-(x(i)-x(j))^2 / (l^2));
        end
    end

    [Q,G,~] = svd(A);
    G = mi*G;

    
        
elseif c == 5 %Example Matérn(1/2)

    alpha = 1; ni = 1/2;
    mi = 10^2;
    kernel = @(x,y) sqrt(pi)*((alpha*norm(x-y))^(ni)*besselk(ni,alpha*norm(x-y)))/(2^(ni-1)*alpha^(2*ni)*gamma(ni+0.5));
    data_matrix = randn(1,n);
    
    for row = 1:n

        for column = 1:n

            if row == column
            
                A(row,column) = sqrt(pi)*gamma(ni)/(gamma(ni + 1/2)*alpha^(2*ni));

            else

                A(row,column) = kernel(data_matrix(:,row),data_matrix(:,column));

            end

        end
    
    end

    [Q,G,~] = svd(A);
    G=mi*G;
        
    
elseif c == 6 %Example Matérn(3/2)

    alpha = 1; ni = 3/2;
    mi = 10^4;
    kernel = @(x,y) sqrt(pi)*((alpha*norm(x-y))^(ni)*besselk(ni,alpha*norm(x-y)))/(2^(ni-1)*alpha^(2*ni)*gamma(ni+0.5));
    data_matrix = randn(1,n);
    
    for row = 1:n

        for column = 1:n

            if row == column
            
                A(row,column) = sqrt(pi)*gamma(ni)/(gamma(ni + 1/2)*alpha^(2*ni));

            else

                A(row,column) = kernel(data_matrix(:,row),data_matrix(:,column));

            end

        end
    
    end

    [Q,G,~] = svd(A);
    G=mi*G;
        
end

trG = sum(log(diag(G+eye(n,n))),"all");