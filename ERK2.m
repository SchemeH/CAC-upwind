clear
xL=-0.5;xR=0.5; % computational domain
epsilon=0.01; % Interface width
kappa=1;  % Stabilization parameter

MBP=[];  % Store extreme values

% Initialization
UExact=@(x,y,t) 0.9*sin(100*pi*x).*sin(100*pi*y);
VxExact=@(x,y,t) -exp(-t).*cos(2*pi*y); % Velocity component in the x-direction
VyExact=@(x,y,t) exp(-t).*sin(2*pi*x); % Velocity component in the y-direction

TT=50; %Total time
tau=0.1; % Time step size
N=64; % Spatial division
h=(xR-xL)/N; % Spatial step size
  x=xL:h:xR-h;
 [X, Y] = meshgrid(x, x);


  Un=reshape(  UExact(X,Y,0) ,[N.^2,1] ); % Initial value
  [M,K]=GenerateMatrix_MK(N,h); % M is the stiffness matrix, and K is the discrete matrix of the central difference of the Laplace operator

% ERK2 
  Ntau=TT/tau; 
 for tnode=1:Ntau
     t=tnode*tau
     Unn=Un;

  Vx=reshape(VxExact(X,Y,t-tau)+1e-16 ,[N.^2,1] );
  Vy=reshape(VyExact(X,Y,t-tau)+1e-16 ,[N.^2,1] );
  [Ex,Ey]=GenerateMatrix_E(N,h,Vx,Vy);
  L=epsilon^2*(1-Un.^2).*K+abs(Vx).*Ex+abs(Vy).*Ey-kappa*M;
  F1=(1-Un.^2).*(Un-Un.^3)+kappa*Un;
  Un=phipm(tau,L,[Un,F1],eps,true);
  % expmv(L, Unn, 1/2*tau)+1/2*tau*(  expmv(L, (1/2*tau*L )\F,1/2*tau )-(M/(1/2*tau*L ))*F  ); 

  Vx=reshape(VxExact(X,Y,t-tau/2)+1e-16 ,[N.^2,1] );
  Vy=reshape(VyExact(X,Y,t-tau/2)+1e-16 ,[N.^2,1] );
  [Ex,Ey]=GenerateMatrix_E(N,h,Vx,Vy);
  L=epsilon^2.*K.*(1-Un.^2)+abs(Vx).*Ex+abs(Vy).*Ey-kappa*M;
  F=(1-Un.^2).*(Un-Un.^3)+kappa*Un;
  % Un=expmv(L, Unn, tau)+tau*(  expmv(L, (tau*L )\F,tau)-(M/(tau*L ))*F  ); 
  Un = phipm(tau,L,[Un,F1,1/tau*(F - F1)],eps,true);
  MBP=[MBP,max(max(abs(Un)))];


 end


figure
surf(X,Y,reshape(Un,[N,N]))
axis([xL,xR,xL,xR])
shading interp
colormap jet
view(0,90)
grid off
box on


function [M,K]=GenerateMatrix_MK(N,h)

I= ones(N,1);
Ke=spdiags([I -2*I I], -1:1, N, N);
Ke(N,1)=1; 
Ke(1,N)=1;
Ke=Ke/h^2;
Me=speye(N,N);

K=kron(Ke,Me)+kron(Me,Ke);
M=kron(Me,Me);
end


function [Ex,Ey]=GenerateMatrix_E(N,h,Vx,Vy)

I= ones(N*N,1);
Ey=(spdiags([0*I  -2*I  0*I], [-1,0,1], N*N, N*N))+diag(1+sign(Vy(2:end)),-1)+diag(1-sign(Vy(2:end)),1);

 i =(1:N:N^2-N+1)';
 j = i+N-1;
 k=i(2:end)-1;
 Ey(sub2ind(size(Ey), i, j))=1+sign(Vy(i));
 Ey(sub2ind(size(Ey), i(2:end), k))=0;

 i=(N:N:N^2)';
 j=i-(N-1);
 k=i(1:end-1)+1;
 Ey(sub2ind(size(Ey), i, j))=1-sign(Vy(i));
 Ey(sub2ind(size(Ey), i(1:end-1), k))=0;
 Ey=Ey/(2*h);
 Ey=sparse(Ey);

 
Ex=-2*speye(N^2);
i=(N+1:N^2-N)';
j=i+N;
k=i-N;
Ex(sub2ind(size(Ex), i, j))=1-sign(Vx(i));
Ex( sub2ind(size(Ex), i, k) )=1+sign(Vx(i));
i=(N^2-N+1:N^2)';
j=(1:N)';
k=i-N;
Ex(sub2ind(size(Ex), i, j))=1-sign(Vx(i));
Ex( sub2ind(size(Ex), i, k) )=1+sign(Vx(i));
i=(1:N)';
j=N^2-N+i;
k=i+N;
Ex(sub2ind(size(Ex), i, j))=1+sign(Vx(i));
Ex( sub2ind(size(Ex), i, k) )=1-sign(Vx(i));
Ex=Ex/(2*h);
end






