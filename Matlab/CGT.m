%% CG
n = 5;
A = hilb(n);
b = rand(n,1);
nm = 100000;%max itertaion

th=1e-5;
x = zeros(n,1);
r = b-A*x;
p = r;
rdr = dot(r,r);

for it = 1:nm
    xprev = x;
    
    alpha = rdr/dot(p,A*p);
    x = x+alpha*p;
    rdrp = rdr;
    r = r - alpha*A*p;
    rdr = dot(r,r);
    beta = rdr/rdrp;
    p = r + beta*p;
    
    inc = norm(xprev-x,inf);
    if inc<th
        break;
    end
end

%% PCG
n = 5;
A = hilb(n);
b = rand(n,1);
nm = 100000;%max itertaion

th = 1e-5;
L = -tril(A,-1);
U = -triu(A,1);
D = diag(diag(A));
omegas = 1.0:0.05:1.9;


%%%%
omega = omegas(3);
S = 1/sqrt(omega*(2-omega))*(D-omega*L)/sqrtm(D);
M = S*S';

x = zeros(n,1);
r = b-A*x;
z = M\r;
p  = z;
zdr = dot(z,r);
for it = 1:nm
    xprev = x;
    alpha = zdr/dot(p,A*p);
    x = x+alpha*p;
    r = r-alpha*A*p;
    z = M\r;
    zdrp = zdr;
    zdr = dot(z,r);
    beta = zdr/zdrp;
    p = z+beta*p;
    
    inc = norm(xprev-x,inf);
    if inc<th
        break;
    end
end
%%%%
