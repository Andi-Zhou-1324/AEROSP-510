clc
clear
close all
%%
coords = [0,0;
          1,0;
          2,2;
          0,1];
element(1,-1,coords)

intpoint = [-1/sqrt(3), -1./sqrt(3);
            -1/sqrt(3), 1/sqrt(3);
            1/sqrt(3), -1/sqrt(3);
            1/sqrt(3), 1/sqrt(3)];
out = 0;
for i = 1:4
    [Ni, J, B] = element(intpoint(i,1), intpoint(i,2), coords);
    temp = Ni(1).*Ni(2).*det(J);
    out = out + temp;
end

%Line Integral
[Ni, J, B] = element(1, 0, coords);
out = Ni(2).*det(J);

out = 2*out;

%Q3e
intpoint = [0.5,0];
[Ni, J, B] = element(0.5, 0, coords);

2.*0.375.*sqrt((0.75-0.375)^2 + (0.25 - 0.875)^2)



%Q2;


%% Functions

function [Ni, J, B] = element(xi, eta, coords) %hw6, p1
    Ni = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    N = [Ni(1) 0 Ni(2) 0 Ni(3) 0 Ni(4) 0
         0 Ni(1) 0 Ni(2) 0 Ni(3) 0 Ni(4) ];

    dNdxi = 0.25*[ eta-1 1-eta 1+eta -1-eta ; xi-1 -1-xi xi+1 1-xi ];
    J = dNdxi*coords;
    dN = J \ dNdxi;
    B = [dN(1, 1) 0 dN(1, 2) 0 dN(1, 3) 0 dN(1, 4) 0
        0 dN(2, 1) 0 dN(2, 2) 0 dN(2, 3) 0 dN(2, 4)
        dN(2, 1) dN(1, 1) dN(2, 2) dN(1, 2) dN(2, 3) dN(1, 3) dN(2, 4) dN(1, 4)
    ];
end

