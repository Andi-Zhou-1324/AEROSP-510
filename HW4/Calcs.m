clc
clear
%%
syms x u_2
x_1 = 0;
x_2 = 0.5;
x_3 = 1;

N_1_1 = (x_2 - x)./(x_2 - x_1);
N_1_2 = (x - x_1)./(x_2 - x_1);

N_2_1 = (x_3 - x)./(x_3 - x_2);
N_2_2 = (x - x_2)./(x_3 - x_2);

diff(N_1_2,x)

eq1 = (int(diff(N_1_2,x)^2,x,x_1,x_2) + int(diff(N_2_1,x)^2,x,x_2,x_3));
eq2 = int(x*N_1_2,x,x_1,x_2) + int(x*N_2_1,x,x_2,x_3)

eq2./eq1