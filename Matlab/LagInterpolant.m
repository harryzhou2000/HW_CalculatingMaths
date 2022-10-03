function yi = LagInterpolant(x,y,xi)
sz1 = size(xi,1);
yi = zeros(numel(xi),1);
xi = reshape(xi,[],1);
x = reshape(x,1,[]);
for i=1:numel(x)
    up = xi - x;
    up(:,i) = 1;
    upp = prod(up,2);
    lo = x-x(i);
    lo(i) = 1;
    lop = prod(lo);
    yi = yi+y(i)/lop*upp;
end
yi = reshape(yi,sz1,[]);