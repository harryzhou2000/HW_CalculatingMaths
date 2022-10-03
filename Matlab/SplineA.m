
%xi = sym([-3, -2, 1, 4 ]);
%fi = sym([2,0,3,1]);
xi = ([-3, -2, 1, 4 ]);
fi = ([2,0,3,1]);

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

A = diag(lai(1:end-1),1)+diag(mui(2:end),-1)+diag(2*ones(1,n-2));
b = di;
b(1) = di(1)- mui(1)*m0;
b(end) = di(end) - lai(end)*mn;

m = [m0;A\b;mn];

OUTS = [xi(thisx),xi(thisx+1),fi(thisx),fi(thisx+1),m(thisx),m(thisx+1)];
hs = cell(1,n-1);
hsf = cell(1,n-1);
hsd = cell(1,n-1);
hsdd = cell(1,n-1);
samps =cell(1,n-1);
for in = 1:n-1
    [hs{in},hsf{in},hsd{in},hsdd{in}] = Hermit3(OUTS(in,:));
    samps{in} = double(linspace(xi(in),xi(in+1),100));
end


%%
hold on;
for in = 1:n-1
    plot(samps{in},hsdd{in}(samps{in}));
end
hold off;




