clc
clear
close all
%%
n_1 = [1,-0.2];
n_2 = [4,0];
n_3 = [1,0.2];


nu = 0.3;
E  = 70E9;

D = E./(1-nu^2).*[1, nu, 0;
                nu, 1, 0;
                0, 0, (1-nu)./2];
t = 10/1E3;

syms xi eta 

N_1 = xi;
N_2 = eta; 
N_3 = 1-xi-eta;

N = [N_1, 0, N_2, 0, N_3, 0;
     0, N_1, 0, N_2, 0, N_3];

rho = 2700;
w = 1000;
%% Computation B, K
x = [n_1(1);n_2(1);n_3(1)];
y = [n_1(2);n_2(2);n_3(2)];

A = (1/2).*(x(1).*(y(2)-y(3)) + x(2).*(y(3)-y(1)) + x(3).*(y(1)-y(2)));


x_13 = x(1) - x(3);
y_23 = y(2) - y(3);

x_23 = x(2) - x(3);
y_13 = y(1) - y(3);

B = [y(2)-y(3), 0, y(3)-y(1), 0, y(1)-y(2), 0;
     0, x(3)-x(2), 0, x(1)-x(3), 0, x(2)-x(1);
     x(3)-x(2), y(2)-y(3), x(1)-x(3), y(3)-y(1), x(2)-x(1), y(1)-y(2)];
B = 1./(x_13.*y_23 - x_23.*y_13).*B;

K = B'*D*B.*t.*A;

%% Change the following for force calculations
F_b = [rho.*w.^2.*(N_1.*x(1) + N_2.*x(2) + N_3.*x(3));
       0];

F = double(int(int(N'*F_b.*(2.*A).*t,eta,0,1-xi),xi,0,1))

K(3:4,3:4)\F(3:4)