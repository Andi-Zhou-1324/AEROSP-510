clc
clear
close all
%%

E = 70E9;
v = 0.3;
quad_coord = [0,0;
              2,0;
              2,1;
              0,2];

ref_coord = [-1/sqrt(3), -1/sqrt(3);
              1/sqrt(3), -1/sqrt(3);
              1/sqrt(3),  1/sqrt(3);
             -1/sqrt(3),  1/sqrt(3)];

surf_coord = [1, 1/sqrt(3);
              1, -1/sqrt(3)];

f = [0;1E6]; %Pa
T = [1E6;0]; %Pa

[N,J,dN] = element(ref_coord(2,1),ref_coord(2,2),quad_coord);

f = localbody(ref_coord,quad_coord)
%% Functions Declared
function [N, J, dN] = element(xi, eta, coords)
    N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dNdxi = 0.25*[ eta-1 1-eta 1+eta -1-eta ; xi-1 -1-xi xi+1 1-xi ];
    J = dNdxi*coords;
    dN = J \ dNdxi;
end

function f = localbody(intpts,coords)
    f = zeros(8,1);
    t = 1;
    for i = 1:size(intpts,1)
        xi = intpts(i,1);
        eta= intpts(i,2);

        N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
        [~, J, ~] = element(xi, eta, coords);
        
        N = [N(1), 0, N(2), 0, N(3), 0, N(4), 0;
            0, N(1), 0, N(2), 0, N(3), 0, N(4)];
        
        f = f + N'*[0;1e6].*det(J).*t;
    end
end

function f = localtraction(tr_node,coords)
    f = zeros(8,1);
    if ~isempty(tr_node) %if nodes have traction
        t = 1;
        intpts = [1 1/sqrt(3);1 -1/sqrt(3)];
        for i = 1:size(intpts,1)
            [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
            detJstar = sqrt(J(2,1)^2 + J(2,2)^2);
            f = f + N'*[1e6;0]*detJstar*t;
        end
    end
end

function k = localstiffnessmat(intpts,coords)
    k = zeros(8,8);
    t = 1;
    E = 70e9;nu = 0.3;
    D = E/(1+nu)/(1-2*nu) * [ 1-nu nu 0 ; nu 1-nu 0 ; 0 0 1/2-nu ];
    for i = 1:size(intpts,1)
        [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
        k = k + B'*D*B*det(J)*t;
    end
end