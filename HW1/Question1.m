clc
clear
close all
%% Question 1
syms k E A L alpha d_T delta_A delta_B P

X = [k + 2*k, -2.*k;
     -2.*k, 3*k + 2*k + E*A/L];

B = [delta_A;
     delta_B];

% F = [0; -E*A*alpha*d_T];

eq1 = X(1,1).*B(1) + X(1,2).*B(2) == 0;
% 
eq2 = X(2,1).*B(1) + X(2,2).*B(2) == -E*A*alpha*d_T + P;

[delta_A,delta_B] = solve(eq1,eq2,delta_A,delta_B)