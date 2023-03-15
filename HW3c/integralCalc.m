clc
clear
%%
syms rho omega w t x x_1 x_2 L

N_1 = (x_2 - x)./(x_2-x_1);
N_2 = (x - x_1)./(x_2-x_1);


eq1 = rho*w^2*x*N_1*t^2;
eq2 = rho*w^2*x*N_2*t^2;


eq1 = (simplify(int(eq1,x,x_1,x_2)));
eq2 = (simplify(int(eq2,x,x_1,x_2)));

eq1_int = subs(eq1,[x_1,x_2],[L/2,L])
eq2_int = subs(eq2,[x_1,x_2],[L/2,L])