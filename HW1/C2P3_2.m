%%
NM = 20;
res = nan(1,NM);
errs = nan(1,NM);
xs = cell(1,NM);
es = cell(1,NM);
for i = 1:NM
   H = hilb(i);
   x = (1:i)';
   b = H*x;
   xbar = H\b;
   res(i) = norm(H*xbar-b,inf);
   xs{i} = xbar;
   es{i} = (xbar-x)./x;
   errs(i) = norm(x-xbar,inf);
end
subplot(1,2,1);
semilogy(res);
xlabel('n');
ylabel('inf norm of residual');
subplot(1,2,2);
semilogy(errs);
xlabel('n');
ylabel('inf norm of error');