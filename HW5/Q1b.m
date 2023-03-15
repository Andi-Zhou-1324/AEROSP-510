clc
clear
close all
%% Question 1B
Q=4000;
k=0.8;
h=20;
T_amb=30;
xL=0;
xR=12.5e-2;

nElem=10;
nNodes=2*nElem+1;
globalA=zeros(nNodes,nNodes);
globalb=zeros(nNodes,1);
dx=(xR-xL)/nElem;
x1=xL;
x3=xL+dx;
x2=(x1+x3)/2;

syms x
N1=((x-x2)./(x1-x2)).*((x-x3)./(x1-x3));
N2=((x-x1)./(x2-x1)).*((x-x3)./(x2-x3));
N3=((x-x1)./(x3-x1)).*((x-x2)./(x3-x2));
N1p=diff(N1,x);
N2p=diff(N2,x);
N3p=diff(N3,x);
f=Q;

for i=1:1:(nElem)
    A11=int(N1p.^2,x,x1,x3);
    A12=int(N1p.*N2p,x,x1,x3);
    A13=int(N1p.*N3p,x,x1,x3);
    A22=int(N2p.^2,x,x1,x3);
    A23=int(N2p.*N3p,x,x1,x3);
    A33=int(N3p.^2,x,x1,x3);

    b11=int(f.*N1,x,x1,x3);
    b21=int(f.*N2,x,x1,x3);
    b31=int(f.*N3,x,x1,x3);

    localA=[A11 A12 A13;
    A12 A22 A23;
    A13 A23 A33];
    
    globalA(2*i-1:2*i+1,2*i-1:2*i+1)=globalA(2*i-1:2*i+1,2*i-1:2*i+1)+localA;
    localb=[b11; b21; b31];
    globalb(2*i-1:2*i+1,1)=globalb(2*i-1:2*i+1,1)+localb;
    x1=x1+dx;
    x2=x2+dx;
    x3=x3+dx;
    disp(localA);
    disp(localb);
    N1=((x-x2)./(x1-x2)).*((x-x3)./(x1-x3));
    N2=((x-x1)./(x2-x1)).*((x-x3)./(x2-x3));
    N3=((x-x1)./(x3-x1)).*((x-x2)./(x3-x2));
    N1p=diff(N1,x);
    N2p=diff(N2,x);
    N3p=diff(N3,x);
    % disp(localb);
end
%BC
A=k*globalA;
b=globalb;
b(end,1)=b(end,1)+T_amb*h;
A(end,end)=A(end,end)+h;
u=A\b;
disp(u);
plot(xL:dx/2:xR,u,'LineWidth',3);
