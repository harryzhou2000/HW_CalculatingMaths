
nm = 1000;
th = 1e-10;
N = 8;
%%

H = hilb(N);
p = diag(min(H,[],2));
%H = p\H;

xa = ones(N,1);
b = H*xa;
L = tril(H);
L = -L.*(ones(N)-eye(N));
U = triu(H);
U = -U.*(ones(N)-eye(N));
D = diag(diag(H));
BJ = eye(N) - D\H;
fJ = D\b;

x = 0*ones(N,1);
for it = 1:nm
    x = BJ*x + fJ;
    res = H*x - b;
    if(norm(res,inf)<1e-10)
        break;
    end
end
x
