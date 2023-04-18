clc
clear
close all
%%
t_vec = [0,90,-30,30];

calculateD(t_vec)



%% Functions Declared

function [C] = calculateC(t) 
    E1 = 145.88E3;
    v12 = 0.263;
    E2 = 13.314E3;
    v23 = 0.47;
    G12 = 4.386E3;
    G23 = 4.528E3;
    E3 = E2;
    G13 = G12;
    v13 = v12;
    
    S_0 = [1/E1, -v12/E1, -v13/E1, 0, 0, 0;
           -v12/E1, 1/E2, -v23/E2, 0, 0, 0;
           -v13/E1, -v23/E2, 1/E3, 0, 0, 0;
           0, 0, 0, 1/G23, 0, 0;
           0, 0, 0, 0, 1/G13, 0;
           0, 0, 0, 0, 0, 1/G12];
    
    l = [cosd(t); cosd(90 + t); cosd(90)];
    m = [cosd(90 -t); cosd(t) ; cosd(90)];
    n = [cosd(90); cosd(90); cosd(0)];
    
    T = [l(1).^2, m(1).^2, n(1).^2, m(1)*n(1), l(1)*n(1), l(1)*m(1); 
        l(2).^2, m(2).^2, n(2).^2, m(2)*n(2), l(2)*n(2), l(2)*m(2);
        l(3).^2, m(3).^2, n(3).^2, m(3)*n(3), l(3)*n(3), l(3)*m(3);
        2*l(2)*l(3), 2*m(2)*m(3), 2*n(2)*n(3), m(2)*n(3)+n(2)*m(3), l(2)*n(3) + n(2)*l(3), l(2)*m(3)+m(2)*l(3);
        2*l(1)*l(3), 2*m(1)*m(3), 2*n(1)*n(3), m(1)*n(3)+n(1)*m(3), l(1)*n(3) + n(1)*l(3), l(1)*m(3)+m(1)*l(3);
        2*l(1)*l(2), 2*m(1)*m(2), 2*n(1)*n(2), m(1)*n(2)+n(1)*m(2), l(1)*n(2) + n(1)*l(2), l(1)*m(2)+m(1)*l(2)];
    
    C_0 = inv(S_0);
    C = T'*C_0*T;
end

function [D] = calculateD(t_vec)
    C = zeros(6);
    for i = 1:size(t_vec,2)
        C = C + (1/size(t_vec,2)).*calculateC(t_vec(i));
    end
    
    S = inv(C);
    
    S(3:5,:) =[];
    S(:,3:5) = [];

    D = inv(S);
end