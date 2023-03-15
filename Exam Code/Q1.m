clc
clear
close all
%%
syms E_1 A_1 L_1 E_2 A_2 L_2 P 

K = [E_1*A_1/L_1 + E_2*A_2/L_2, -E_2*A_2/L_2;
     -E_2*A_2/L_2, E_2*A_2/L_2];
F = [P; -P];

u = K\F;