function yi = SplineNormal(x,y,xi) 
sz1 = size(xi,1);
x = reshape(x,1,[]);
y = reshape(y,1,[]);
xi = reshape(xi,[],1);

this = 1:(numel(x)-1);
h = x(this+1)-x(this); 
thish = 1:(numel(h)-1);
mu = (h(thish+1)./(h(thish)+h(thish+1)));
la = 1-mu;

fd = (y(this+1)-y(this))./h;
d = 3*(mu.*fd(thish) + la.*fd(thish+1));
d = [3*fd(1),d,3*fd(end)];

A = diag([mu,1],-1) + diag([1,la],1) + diag(2*ones(1,numel(x)));
ms = A\(d');


yi = zeros(numel(xi),1);

for i = 1:numel(xi)
    for j = 1:(numel(x)-1)
        if xi(i)>=x(j) && xi(i)<x(j+1)
            alphaA = (1+2*(xi(i)-x(j))/h(j))*((xi(i)-x(j+1))/h(j))^2;
            alphaB = (1-2*(xi(i)-x(j+1))/h(j))*((xi(i)-x(j))/h(j))^2;
            betaA = (xi(i)-x(j))*((xi(i)-x(j+1))/h(j))^2;
            betaB = (xi(i)-x(j+1))*((xi(i)-x(j))/h(j))^2;
            yi(i) = alphaA*y(j) + alphaB*y(j+1) + betaA*ms(j) + betaB*ms(j+1);
            break;
        end
    end
    if(xi(i) == x(end))
        yi(i) = y(end);
    end
end
yi = reshape(yi,sz1,[]);
