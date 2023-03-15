clc
clear
close all
%% Question 3
E = 100E6;
A_1 = 1500*0.0001; %m^2
A_2 = 2500*0.0001; %m^2
syms x

L = 1;
x_1 = 0;
x_2 = L;
x_3 = 2*L;

rho_g = 0.06; %N/cm^3
rho_g = rho_g/1E6; %converts to N/m^3

N_1_1 = (x_2 - x)./L;
N_1_2 = (x - x_1)./L;

N_2_1 = (x_3 - x)./L;
N_2_2 = (x - x_2)./L;

eq1 = -rho_g*N_1_1 * A_1;

eq1 = double(int(eq1,x,0,L));

eq2 = -rho_g*N_1_2 * A_1;

eq2 = double(int(eq2,x,0,L));

eq3 = -rho_g*N_2_1 * A_2;

eq3 = double(int(eq3,x,L,2*L));

eq4 = -rho_g*N_2_2 * A_2;

eq4 = double(int(eq4,x,L,2*L));

F = [eq2 + eq3;eq4];

K = [E*A_1/L + E*A_2/L, -E*A_2/L;
     -E*A_2/L, E*A_2/L];

u = K\F