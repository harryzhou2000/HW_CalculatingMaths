N = 20;
nm = 100000;
th = 1e-6;
h = 1/(N+1);
IBR = @(i,j,nc) (i-1)*nc + j;
N = N-2;
%A = zeros((N+1)^2);
A = sparse(N+1,N+1);
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
colormap jet
imagesc(A);
colorbar;
axis equal;
title(sprintf('N = %d',N+2));
%%
figure;
plot(b);
title(sprintf('N = %d',N+2));
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


