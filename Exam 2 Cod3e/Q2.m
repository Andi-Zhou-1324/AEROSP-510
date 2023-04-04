clc
clear
close all
%% 6 noded Tri

syms xi eta


N_1 = 1-3*xi+2*xi^2-3*eta+4*eta*xi+2*eta^2;
N_2 = -xi + 2*xi^2;
N_3 = -eta + 2*eta;
N_4 = 4*xi-4*xi^2-4*xi*eta;
N_5 = 4*xi*eta;
N_6 = 4*eta-4*xi*eta - 4*eta^2;
N_mat = [N_1, N_2, N_3, N_4, N_5, N_6];

N = [N_1 0 N_2 0 N_3 0 N_4 0 N_5 0 N_6 0
         0 N_1 0 N_2 0 N_3 0 N_4 0 N_5 0 N_6 ];
coords = [0,0;
          2,0;
          0,2;
          1,0;
          1,1;
          0,1];



intpoint = [1/3, 0;
            2/3, 0];

F = zeros(2*6,1);
for i = 1:2
    detJstar = element_tri(intpoint(i,1),intpoint(i,2),coords);

    N = [N_1 0 N_2 0 0 0 N_4 0 0 0 0 0
         0 N_1 0 N_2 0 0 0 N_4 0 0 0 0 ];
    N = subs(N,[xi,eta],[intpoint(i,1),intpoint(i,2)]);
    N = double(N);
    
    T = [0;(coords(1,1)*N_1 + coords(2,1)*N_2 + coords(3,1)*N_3 + coords(4,1)*N_4 + coords(5,1)*N_5 + coords(6,1)*N_6)];
    T = subs(T,[xi,eta],[intpoint(i,1),intpoint(i,2)]);
    T = double(T);
    
    temp = N'*T.*detJstar;
    F = F + temp;

end

function [detJstar] = element_tri(xi_num, eta_num, coords) %hw6, p1
    syms xi eta

    
    N_1 = 1-3*xi+2*xi^2-3*eta+4*eta*xi+2*eta^2;
    N_2 = -xi + 2*xi^2;
    N_3 = -eta + 2*eta;
    N_4 = 4*xi-4*xi^2-4*xi*eta;
    N_5 = 4*xi*eta;
    N_6 = 4*eta-4*xi*eta - 4*eta^2;
    N_mat = [N_1, N_2, N_3, N_4, N_5, N_6];

   
    
    dNdxi = [diff(N_mat,xi); diff(N_mat,eta)];
    dNdxi = subs(dNdxi,[xi,eta],[xi_num,eta_num]);
    dNdxi = double(dNdxi);
    
    J = dNdxi*coords;
    detJstar = det(J);
%     detJstar = sqrt(J(1,1)^2 + J(1,2)^2); 
end