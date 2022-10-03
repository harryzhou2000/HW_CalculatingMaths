
xi = sym([-3, -2, 1, 4 ]);
fi = sym([2,0,3,1]);

xi = reshape(xi,[],1);
fi = reshape(fi,[],1);
m0 = -1;
mn = 1;
n = numel(xi);

thisx = 1:(n-1);
thish = 1:(n-2);
hi = xi(thisx+1)-xi(thisx);
mui = hi(thish+1)./(hi(thish) + hi(thish+1));
lai = hi(thish)./(hi(thish) + hi(thish+1));
fp = (fi(thisx+1)-fi(thisx))./hi;
di = 3*(fp(thish).*mui + lai.*fp(thish+1));

A = diag([1;lai],1)+diag([mui;1],-1)+diag(2*ones(1,n));
b = [3*fp(1);di;3*fp(end)];
% b(1) = di(1) - mui(1)*m0;
% b(end) = di(end) - lai(end)*mn;

m = A\b;

OUTS = [xi(thisx),xi(thisx+1),fi(thisx),fi(thisx+1),m(thisx),m(thisx+1)];
hs = cell(1,n-1);
hsf = cell(1,n-1);
hsd = cell(1,n-1);
hsdd = cell(1,n-1);
for in = 1:n-1
    [hs{in},hsf{in},hsd{in},hsdd{in}] = Hermit3(OUTS(in,:));
end


%%
hold on;
for in = 1:n-1
    plot(samps{in},hsf{in}(samps{in}));
end
hold off;


