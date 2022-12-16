clear; clf;
N = 21; % partition number
nm = 100000;
th = 1e-6;
h = 1/(N+1);
IBR = @(i,j,nc) (i-1)*nc + j;
N = N-2;
%A = zeros((N+1)^2);
A = sparse((N+1)^2,(N+1)^2);
b = zeros((N+1)^2,1);
u = b;
ua = b;
for i = 1:N+1
    for j = 1:N+1
        x = (i)*h;
        y = (j)*h;
        I = IBR(i,j,N+1);
        A(I,I) = 4;
        b(I) = h^2 * 2*pi^2*sin(pi*x)*sin(pi*y);
        ua(I) = sin(pi*x)*sin(pi*y);
        if i>1
           J = IBR(i-1,j,N+1);
           A(I,J) = -1;
        end
        if i<=N
           J = IBR(i+1,j,N+1);
           A(I,J) = -1;
        end
        if j>1
           J = IBR(i,j-1,N+1);
           A(I,J) = -1;
        end
        if j<=N
           J = IBR(i,j+1,N+1);
           A(I,J) = -1;
        end
    end
end
% colormap jet
colormap(linspace(1,0,1001)'*ones(1,3));
imagesc(full(abs(A)));
colorbar;
axis equal;
title(sprintf('N = %d',N+1));
set(gca,'FontName','Times New Roman');
%%
%figure;
plot(b);
title(sprintf('N = %d',N+1));
set(gca,'FontName','Times New Roman');
%% plot solution
%x = A\b;
colormap jet
us = A\b;
us = reshape(us,N+1,N+1);
[xg,yg] = meshgrid(linspace(h,1-h,N+1),linspace(h,1-h,N+1));
errmax = max(abs(us - sin(pi*xg).*sin(pi*yg)),[],'all')

surf(xg,yg,us,'LineStyle','none');
hold on;
% surf(xg,yg,sin(pi*xg).*sin(pi*yg) -0.0);
hold off;
colorbar;
view(3);
xlabel('x');ylabel('y');zlabel('u');
set(gca,'FontName','Times New Roman');
%% plot err
surf(xg,yg,abs(us-sin(pi*xg).*sin(pi*yg)),'LineStyle','none');
colorbar;
view(3);
xlabel('x');ylabel('y');zlabel('err');
set(gca,'FontName','Times New Roman');
%% SOR
L = -tril(A,-1);
U = -triu(A,1);
D = diag(diag(A));
omegas = 1.0:0.05:1.9;
errs = cell(numel(omegas),1);
names = cell(numel(omegas),1);
errss = omegas*nan;
itts = omegas*nan;
i = 0;

for omega = omegas
    i = i+1;
    errs{i} = nan(nm,1);
    Bsor = (D-omega*L)\((1-omega)*D+omega*U);
    fsor = omega*((D-omega*L)\b);
    
    x = zeros((N+1)^2,1);
    for it = 1:nm
        xprev = x;
        x = Bsor*x + fsor;
        err = x-ua;
        errn = norm(err,inf);
        resn = norm(A*x-b,inf);
        inc = norm(xprev-x,inf);
        %errs{i}(it) = norm(err,inf);
        if inc<th
            break;
        end
        fprintf('err = %f\n',errn);
    end
    itts(i) = it;
    names{i} = sprintf('SOR \\omega = %f',omega);
    errss(i) = norm(err,inf);
end

t = tiledlayout(1,2,'TileSpacing','Compact');
title(t,sprintf('N = %d',N+1),'FontName','Times New Roman');

nexttile(1);
plot(omegas,errss,'*');
xlabel('SOR \omega');
ylabel('norm_{\infty}(err)');
set(gca,'FontName','Times New Roman');
grid on;
nexttile(2);
plot(omegas,itts,'*');
xlabel('SOR \omega');
ylabel('iterations');
set(gca,'FontName','Times New Roman');
grid on;


a = gcf;
a.Units = 'pixels';
a.Position = [100,100,1000,400];


%% CG
x = zeros((N+1)^2,1);
r = b-A*x;
p = r;
rdr = dot(r,r);

for it = 1:nm
   xprev = x;
   
   alpha = rdr/dot(p,A*p);
   x = x+alpha*p;
   rdrp = rdr;
   r = r - alpha*A*p;
   rdr = dot(r,r);
   beta = rdr/rdrp;
   p = r + beta*p;
    
   inc = norm(xprev-x,inf);
   if inc<th
      break; 
   end
end
err = norm(x-ua,inf)


%% PCG

L = -tril(A,-1);
U = -triu(A,1);
D = diag(diag(A));
omegas = 1.0:0.05:1.9;
errs = cell(numel(omegas),1);
names = cell(numel(omegas),1);
errss = omegas*nan;
itts = omegas*nan;
i = 0;

for omega = omegas
    i = i+1;
   S = 1/sqrt(omega*(2-omega))*(D-omega*L)/sqrtm(D);
   M = S*S';
    
    x = zeros((N+1)^2,1);
    r = b-A*x;
    z = M\r;
    p  = z;
    zdr = dot(z,r);
    for it = 1:nm
        xprev = x;
        alpha = zdr/dot(p,A*p);
        x = x+alpha*p;
        r = r-alpha*A*p;
        z = M\r;
        zdrp = zdr;
        zdr = dot(z,r);
        beta = zdr/zdrp;
        p = z+beta*p;
        
        inc = norm(xprev-x,inf);
        if inc<th
           break; 
        end
    end
    itts(i) = it;
    errss(i) = norm(x-ua,inf);
end

t = tiledlayout(1,2,'TileSpacing','Compact');
title(t,sprintf('N = %d',N+1),'FontName','Times New Roman');

nexttile(1);
plot(omegas,errss,'*');
xlabel('PCG SSOR preconditioner \omega');
ylabel('norm_{\infty}(err)');
set(gca,'FontName','Times New Roman');
grid on;

nexttile(2);
plot(omegas,itts,'*');
xlabel('PCG SSOR preconditioner \omega');
ylabel('iterations');
yticks(min(itts):max(itts));
set(gca,'FontName','Times New Roman');
grid on;


a = gcf;
a.Units = 'pixels';
a.Position = [100,100,1000,400];






    