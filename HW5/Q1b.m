%AE 510 Class code for live lecture 
%Author: Your instructor
 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex')
clc
clear
close all
%%

Len = 12.5/100;%length of the bar
k = 0.8; %W/m - C
Q = 4000; %W/m^3

nelem =10;%number of elements (keep it an even number for this problem)
h = 20; %W/m^3*C
T_amb = 30; %Celsius

%%%%%%%%%%%%PREPROCESSING%%%%%%%%%%%
%coordinate matrix [x,y] for each node
co = [0 : Len/(2*nelem): Len]';

%element-node connectivity matrix, length, area, modulus
e = [];
for i = 1:2:2*nelem
    temp = [i,i+1,i+2];
    e = [e;temp];
end
 
Nel = size(e,1);%number of elements
Nnodes = size(co,1); %number of nodes
nne = 3; %number of nodes per element
dof = 1; %degree of freedom per node

%%%%%%%%%%%%PREPROCESSING END%%%%%%%%%%%

 
%%%Generic block: Initializes global stiffness matrix 'K' and force vector 'F'
K = zeros(Nnodes*dof,Nnodes*dof);
F = zeros(Nnodes*dof,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%Assemble Global system - generic FE code 
for A = 1:Nel
    
    syms x
    x_co = co(e(A,:),:);
    
    N_1 = ((x - x_co(2))*(x - x_co(3)))./((x_co(1) - x_co(2))*(x_co(1) - x_co(3)));
    N_2 = ((x - x_co(1))*(x - x_co(3)))./((x_co(2) - x_co(1))*(x_co(2) - x_co(3)));
    N_3 = ((x - x_co(1))*(x - x_co(2)))./((x_co(3) - x_co(1))*(x_co(3) - x_co(2)));
    
    dN_1 = diff(N_1,x);
    dN_2 = diff(N_2,x);
    dN_3 = diff(N_3,x);
    localstiffness = [dN_1.^2, dN_1*dN_2, dN_1*dN_3;
                      dN_1*dN_2, dN_2.^2, dN_2*dN_3;
                      dN_1*dN_3, dN_2*dN_3,dN_3.^2];

    localstiffness = double(int(k.*localstiffness,x,x_co(1),x_co(3)));

    localforce = [N_1;N_2;N_3];
    localforce = double(int(localforce.*Q,x,x_co(1),x_co(3)));

    %DONT TOUCH BELOW BLOCK!! Assembles the global stiffness matrix, Generic block which works for any element    
 
    for B = 1: nne
        for i = 1: dof
            nK1 = (e(A, B)-1)*dof+i;
            nKe1 = (B-1)*dof+i;
            F(nK1) = F(nK1) + localforce(nKe1);
            for C = 1: nne
                for j = 1: dof
                    nK2 = (e(A, C)-1)*dof+j;
                    nKe2 = (C-1)*dof+j;
                    K(nK1, nK2) = K(nK1, nK2) + localstiffness(nKe1, nKe2);
                end
            end
        end
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


%%%%%%%%%%%%%%%%%%%BOUNDARY CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%
 
%external forces

F(end) = F(end) + h*T_amb; %given x- component of force in node 2 
K(end,end) = K(end,end) + h;

%Apply displacement BC by eliminating rows and columns of nodes 3-4 (corresponding to
%degrees of freedom 5 to 8) - alternative (and more generic method) is the penalty approach, or
%static condensation approach - see later class notes

% deletedofs = [1];%first nodes
% K(deletedofs,:) = [];
% K(:,deletedofs) = [];
% F(deletedofs,:) = [];

%%%%%%%%%%%%%%%%%%%BOUNDARY CONDITIONS END%%%%%%%%%%%%%%%%%%%%%%%

 
%solve for displacement unknowns (uk)
uk = K\F;

%expand u to include deleted displacement bcs
% u = ones(Nnodes*dof,1);
% u(deletedofs) = 0;
% I = find(u == 1);
% u(I) = uk;


figure()
hold on
plot(co,uk,'LineWidth',1.5)
title('Temperature Distribution along the Plate')
xlabel('x (m)')
ylabel('T(Celsius)')
xlim([0,Len])
grid on 