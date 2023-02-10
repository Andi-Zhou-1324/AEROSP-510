clc
clear
close all
%% Question 1
t = 0.01; %cm
E = 70E9;
A = t^2;
rho = 2700; %kg/m^3
L = 1; %m
omega = 1000; %rad/s

syms b c d x

F_c = rho*A*omega^2*x;
F_point = 0.1*rho*A*omega^2*x;

u = b*x + c*x^2 + d*x^3;

term_1 = 1/2*E*A*diff(u,x)^2;
term_1_int = simplify(int(term_1,x,0,L));

term_2 = subs(F_point,x,L)*subs(u,x,L);

term_3 = u*F_c;
term_3_int = int(term_3,x,0,L);

Pi = term_1_int - term_2 - term_3_int;

eq1 = diff(Pi,b);
eq2 = diff(Pi,c);
eq3 = diff(Pi,d);

[b,c,d] = solve(eq1,eq2,eq3,b,c,d);

b = double(b);
c = double(c);
d = double(d);

%Plotting:

x = linspace(0,L,1000);

u = 0 + b.*x + c.*x.^2 + d.*x.^3;

figure
hold on
plot(x,u,'LineWidth',1.5)
xlabel('$x (m)$','Interpreter','latex')
ylabel('$u(m)$','Interpreter','latex')
title('Rayleigh - Ritz Method','Interpreter','latex')
grid on

%%  Potential Energy Minimization 
syms d_1 d_2

m_half = 1/2*rho*A*L;
m_full = 1*rho*A*L;

k = 2*E*A/L;

Pi = 1/2*k*d_1^2 + (1/2)*k*(d_1+d_2)^2 - (m_half * omega^2*(d_1) + m_full*omega^2*(d_1 + d_2)) - (0.1*rho*L*t^2*omega^2*L)*(L+d_1+d_2);

dPi_d1 = diff(Pi,d_1);
dPi_d2 = diff(Pi,d_2);

[d_1_solve,d_2_solve] = solve(dPi_d1,dPi_d2,d_1,d_2);


d_1_solve = double(d_1_solve);
d_2_solve = double(d_2_solve);

plot([0,0.5,1],[0,d_1_solve,d_2_solve],'-o','LineWidth',1.5)
legend('Rayleigh - Ritz', 'PE Minimization','Location','northwest')