clc
clear
close all
%% 
E = 70E9;
w = 1000;
L = 1;
A = 0.01^2;
rho = 2700;
M = 0.1.*rho*L*A

M*w^2*L./(E*A) + rho*w^2./E*(1/3)