% clear; close all;
f1 = @(x) 1./(1+25*x.^2);
f2 = @func2;
f = f1;
f = f2;

t = tiledlayout(2,2,'TileSpacing','compact');




%% lags
nexttile(1);
xs = -1:0.1:1;
ys = f(xs);
xi = -1:0.01:1;
plot(xi,f(xi),'k',xi,LagInterpolant(xs,ys,xi),'-.b',xs,ys,'ok');
legend('Precise','Interp','Nodes','Location','best');
title('Newtonian','FontName','Times New Roman');
set(get(gca,'legend'),'FontName','Times New Roman');
ylim([min(f(xi))-0.2,1.5]);
%% lags cheb
nexttile(2);
xs = cos((2*(1:21)-1)/42*pi);
ys = f(xs);
xi = -1:0.01:1;
plot(xi,f(xi),'k',xi,LagInterpolant(xs,ys,xi),'-.b',xs,ys,'ok');
legend('Precise','Interp','Nodes','Location','best');
title('Lagrangian (Chebyshev Nodes)','FontName','Times New Roman');
set(get(gca,'legend'),'FontName','Times New Roman');
ylim([min(f(xi))-0.2,1.5]);
%% spline N

nexttile(3);
xs = -1:0.1:1;
ys = f(xs);
xi = -1:0.01:1;
plot(xi,f(xi),'k',xi,SplineNormal(xs,ys,xi),'-.b',xs,ys,'ok');
legend('Precise','Interp','Nodes','Location','best');
title('Cubic Spline','FontName','Times New Roman');
set(get(gca,'legend'),'FontName','Times New Roman');
ylim([min(f(xi))-0.2,1.5]);
%% linear


nexttile(4);
xs = -1:0.1:1;
ys = f(xs);
xi = -1:0.01:1;
plot(xi,f(xi),'k',xi,LinearInterpolant(xs,ys,xi),'-.b',xs,ys,'ok');
legend('Precise','Interp','Nodes','Location','best');
title('Piecewise Linear','FontName','Times New Roman');

set(get(gca,'legend'),'FontName','Times New Roman');
ylim([min(f(xi))-0.2,1.5]);


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

%%

