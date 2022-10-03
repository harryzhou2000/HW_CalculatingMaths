
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



omegas = 1.0:0.25:1.75;
errs = cell(numel(omegas),1);
names = cell(numel(omegas),1);
errss = omegas*nan;
i = 0;



for omega = omegas
    i = i+1;
    errs{i} = nan(nm,1);
    Bsor = (D-omega*L)\((1-omega)*D+omega*U);
    fsor = omega*((D-omega*L)\b);
    
    x = 0*ones(N,1);
    for it = 1:nm
        x = Bsor*x + fsor;
        res = H*x - b;
        err = x-xa;
        errs{i}(it) = norm(x-xa,inf);
    end
    if(i>1)
        hold on;
    end
    semilogy(errs{i});
    names{i} = sprintf('\\omega = %f',omega);
    errss(i) = errs{i}(end);
end

legend(names);

title(sprintf('n = %d',N));
xlabel('iter');
ylabel('max error');
hold off;