N = 20;
A = zeros(N,N);
A = A + diag (1:N) * 2;
A = A + diag(ones(N-1,1), 1) * 1.5;
A = A + diag(ones(N-1,1), -1) * 14;
b = A*ones(N,1);

[x, flag, relres, iter, resvec] = gmres(A,b,3,[],100,[],[],zeros(N,1));
resvec/resvec(1)