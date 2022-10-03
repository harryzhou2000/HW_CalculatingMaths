function [out,outf,outd,outdd] = Hermit3(s1)
syms x xk xkp fk fkp fdk fdkp
alphak = (1+2*(x-xk)/(xkp-xk))*((x-xkp)/(xk-xkp))^2;
alphakp = (1+2*(x-xkp)/(xk-xkp))*((x-xk)/(xkp-xk))^2;
betak = (x-xk)*((x-xkp)/(xk-xkp))^2;
betakp = (x-xkp)*((x-xk)/(xkp-xk))^2;
% 
bess = alphak*fk + alphakp*fkp + betak * fdk + betakp *fdkp;
bessdd = diff(bess,x,2);
bessddk = simplify(subs(bessdd,xk));
bessddkp = simplify(subs(bessdd,xkp));



bessact = expand(subs(bess,[xk,xkp,fk,fkp,fdk,fdkp],s1));
out = bessact;
outf = matlabFunction(symfun(bessact,x));
outd = matlabFunction(symfun(diff(bessact,1),x));
outdd = matlabFunction(symfun(diff(bessact,2),x));



