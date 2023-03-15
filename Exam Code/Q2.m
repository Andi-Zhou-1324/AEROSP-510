clc
clear
close all
%%
syms E A_1 A_2 L k P R_1 k_bar

A_avg = (A_1 + A_2)/2;

K_avg = E*A_avg/L;

k_bar = simplify((1/3)*(E*A_avg/L) + (1/3)*(E*A_1/L) + (1/3)*(E*A_2/L));


u_2 = simplify(-2*P./(k_bar + k))