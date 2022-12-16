clear;
nm = 1000;
th = 1e-10;
N = 12;
%%

H = hilb(N);
p = diag(min(H,[],2));
%H = p\H;

xa = ones(N,1);
b = H*xa;
L = -tril(H,-1);
U = -triu(H,1);
D = diag(diag(H));
BJ = eye(N) - D\H;
fJ = D\b;

%%

x = 0*ones(N,1);
for it = 1:nm
    x = BJ*x + fJ;
    res = H*x - b;
    if(norm(res,inf)<1e-10 ||   sum(isnan(x)))
        break;
    end
end
x


%%
omegas = 1.0:0.25:1.75;
errs = cell(numel(omegas),1);
ress = cell(numel(omegas),1);
names = cell(numel(omegas),1);
errss = omegas*nan;
resss = omegas*nan;
i = 0;

figure(111);
clf;
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','none');
title(t,sprintf('n = %d',N),'FontName','Times New Roman');



for omega = omegas
    i = i+1;
    errs{i} = nan(nm,1);
    ress{i} = nan(nm,1);
    opt.LT = true;
    Bsor = linsolve(D-omega*L, ((1-omega)*D+omega*U), opt);
    %     Bsor = (D-omega*L)\((1-omega)*D+omega*U);
    fsor = omega*((D-omega*L)\b);
    
    x = 0*ones(N,1);
    for it = 1:nm
        x = Bsor*x + fsor;
        res = H*x - b;
        err = x-xa;
        errs{i}(it) = norm(err, inf);
        ress{i}(it) = norm(res, inf);
    end
    
    
    
    figure(111);
    nexttile(1);
    hold on;
    semilogy(errs{i});
    errss(i) = errs{i}(end);
    nexttile(2);
    hold on;
    semilogy(ress{i});
    resss(i) = ress{i}(end);
    
    names{i} = sprintf('\\omega = %.2f',omega);
end


figure(111);
nexttile(1);
legend(names);

xlabel('iter');
ylabel('norm_{\infty}(err)');
hold off;
set(gca,'FontName','Times New Roman');
set(gca,'YScale','log');
grid on;

nexttile(2);
legend(names);
xlabel('iter');
ylabel('norm_{\infty}(res)');
hold off;
set(gca,'FontName','Times New Roman');
set(gca,'YScale','log');
grid on;