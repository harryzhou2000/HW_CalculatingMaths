%%
NM = 12;
conds = nan(1,NM);
for i = 2:NM
   H = hilb(i);
   conds(i) = cond(H,inf); 
end


lc = log(conds);
[a,~,~,~,stats] = regress(lc',[(1:NM)',(1:NM)'.^0]);
% reg
ns = linspace(1,NM,1001);
condsEst = exp(polyval(a, ns));

plot(1:NM,conds,'O', ns, condsEst);
set(gca,'YScale','log');
xlabel('n');
ylabel('cond_{\infty}');
grid on; grid minor;
title('Hilbert Matrix cond_{\infty}');
set(gca,'FontName','Times New Roman');
legend('Data','fit');

%%
NM = 20;
res = nan(1,NM);
errs = nan(1,NM);
xs = cell(1,NM);
es = cell(1,NM);
for i = 1:NM
   H = hilb(i);
%    x = (1:i)';
   x = ones(i,1);
   b = H*x;
   opt.SYM = true;
   xbar = linsolve(H,b,opt);
%    R = chol(H);
%    xbar = R\(R'\b);
   res(i) = norm(H*xbar-b,inf);
   xs{i} = xbar;
   es{i} = (xbar-x)./x;
   errs(i) = norm(x-xbar,inf);
end
subplot(1,2,1);
semilogy(res,'o');
xlabel('n');
ylabel('norm_{\infty}(res)');
grid on;
set(gca,'FontName','Times New Roman');
xlim([0,NM]);
subplot(1,2,2);
semilogy(errs,'o');
xlabel('n');
ylabel('norm_{\infty}(err)');
grid on;
set(gca,'FontName','Times New Roman');
xlim([0,NM]);