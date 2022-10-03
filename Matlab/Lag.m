function y = Lag(xi,fxi)

npo = numel(xi);
sz = size(xi);
fxi = reshape(fxi,1,[]);
if sz(1) == 1 
    xxi = repmat(xi',1,npo);
    
else
    xxi = repmat(xi,1,npo);
end
k = ones(npo) - eye(npo);
syms xa;
y = sum(prod((xxi-xa).*k + eye(npo),1)./prod((xxi - xxi').*k + eye(npo),1).*fxi);
y = expand(y);
y = matlabFunction(symfun(y,xa));





