%%
NM = 12;
conds = nan(1,NM);
for i = 1:NM
   H = hilb(i);
   conds(i) = cond(H,inf); 
end
semilogy(conds);
xlabel('n');
ylabel('inf cond');

lc = log(conds);
[a,~,~,~,stats] = regress(lc',[(1:NM)',(1:NM)'.^0]);
