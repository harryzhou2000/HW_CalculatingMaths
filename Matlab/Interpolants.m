f1 = @(x) 1./(1+25*x.^2);
f2 = @func2;
f = f2;

%% lags
subplot(2,2,1);
xs = -1:0.1:1;
ys = f(xs);
xi = -1:0.01:1;
plot(xi,f(xi),'b',xi,LagInterpolant(xs,ys,xi),'--r',xs,ys,'ok');
legend('Original','Interpolation','Nodes','Location','best');
title('Newtonian Interpolation','FontName','Cambria Math');
set(get(gca,'legend'),'FontName','Cambria Math');
%% lags cheb
subplot(2,2,2);
xs = cos((2*(0:20))/42*pi);
ys = f(xs);
xi = -1:0.01:1;
plot(xi,f(xi),'b',xi,LagInterpolant(xs,ys,xi),'--r',xs,ys,'ok');
legend('Original','Interpolation','Nodes','Location','best');
title('Lagrangian Interpolation with Chebyshev Nodes','FontName','Cambria Math');
set(get(gca,'legend'),'FontName','Cambria Math');
%% spline N

subplot(2,2,3);
xs = -1:0.1:1;
ys = f(xs);
xi = -1:0.01:1;
plot(xi,f(xi),'b',xi,SplineNormal(xs,ys,xi),'--r',xs,ys,'ok');
legend('Original','Interpolation','Nodes','Location','best');
title('Cubic Spline Interpolation','FontName','Cambria Math');
set(get(gca,'legend'),'FontName','Cambria Math');
%% linear


subplot(2,2,4);
xs = -1:0.1:1;
ys = f(xs);
xi = -1:0.01:1;
plot(xi,f(xi),'b',xi,LinearInterpolant(xs,ys,xi),'--r',xs,ys,'ok');
legend('Original','Interpolation','Nodes','Location','best');
title('Linear Interpolation','FontName','Cambria Math');

set(get(gca,'legend'),'FontName','Cambria Math');

%%



function y = func2(x)
    yA = sin(pi*x);
    yA(x>=0) = 0;
    yB = cos(pi*x);
    yB(x<0) = 0;
    yB(x>=0.5) = 0;
    %y = zeros(size(x,1),size(x,2));
    y = yA + yB;
end