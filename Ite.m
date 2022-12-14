clc;
% fp = @(x) [1-0.7*cos(x(1)),0.2*sin(x(2));0.7*sin(x(1)),1+0.2*cos(x(2)) ];
% f = @(x)[x(1)-0.7*sin(x(1))-0.2*cos(x(2)),x(2) - 0.7*cos(x(1))+0.2*sin(x(2))]';
f = @(x) [x(1)^2 + x(2)^2-4; x(1)^2-x(2)^2-1]
fp = @(x) [2*x(1), 2*x(2);2*x(1),-2*x(2)]

i = 0
x = [1.6,1.2]'
f(x)
fp(x)
dx = - fp(x)\f(x)

%%
clc
i = i+1
x = x - fp(x)\f(x)
f(x)
fp(x)
dx = - fp(x)\f(x)
