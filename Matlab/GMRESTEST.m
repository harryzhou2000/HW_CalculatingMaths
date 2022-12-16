N = 20;
A = zeros(N,N);
A = A + diag (1:N) * 2;
A = A + diag(ones(N-1,1), 1) * 14;
A = A + diag(ones(N-1,1), -1) * 1.5;
b = A*ones(N,1);
AP = A;
bP = b;

D = diag(diag(A));
LL =- tril(A,-1);
UU =- triu(A,1);

ptype = 5;

switch ptype
    case 1
        L = D;
        R = eye(N);
    case 2
        L = sqrt(D);
        R = sqrt(D);
    case 3
        L = eye(N);
        R = D;
    case 4
        L = (D-LL) * D * (D-UU);
        R = eye(N);
    case 5
        L = eye(N);
        R = (D-LL) * D * (D-UU);
    case 6
        L = (D-LL) * sqrt(D);
        R = sqrt(D) * (D-UU);
    case 7
        L = (D-LL);
        R = D * (D-UU);
end







AP = L\A/R;
bP = L\b;



[x, flag, relres, iter, resvec] = gmres(AP,bP,3,1e-15,40,[],[],zeros(N,1));
log10(resvec/resvec(1))

x = R\x;

err = norm(x-ones(N,1),inf)


%%
N = 20;
A = zeros(N,N);
A = A + diag (1:N) * 2;
A = A + diag(ones(N-1,1), 1) * 14;
A = A + diag(ones(N-1,1), -1) * 1.5;
b = A*ones(N,1);
AP = A;
bP = b;



[x, flag, relres, iter, resvec] = gmres(AP,bP,3,1e-15,100,[],[],zeros(N,1));
log10(resvec/resvec(1))
%%
load west0479;
A = west0479;
Aeigs = eigs(A);
N = size(A,1);
A = A + speye(N) * max(abs(Aeigs))*1e-6;
b = A*ones(N,1);
[L,U] = ilu(A);

% AP = U\(L\A);
% bP = U\(L\b);
% AP = (L\A)/U;
% bP = (L\b);
AP = (A/U)/L;
bP =b;

[x, flag, relres, iter, resvec] = gmres(AP,bP,3,1e-15,10,[],[],zeros(N,1));
log10(resvec/resvec(1))