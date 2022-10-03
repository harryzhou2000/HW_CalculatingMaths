function yi = LinearInterpolant(x,y,xi)

this = 1:numel(x)-1;
hs = x(this+1)-x(this);
yi = zeros(numel(xi),1);
for i = 1:numel(xi)
    for j = 1:numel(x)-1
        if xi(i)>=x(j) && xi(i)<x(j+1)
        
            yi(i) = (xi(i)-x(j))/hs(j) * y(j+1) + (x(j+1)-xi(i))/hs(j) * y(j);
            break;
        end
    end
    if(xi(i) == x(end))
       yi(i) = y(end); 
    end
end
yi = reshape(yi,size(xi,1),[]);