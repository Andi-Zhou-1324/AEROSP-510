clc
clear
close all


N = 10000;
[co,e] = buildMesh(N);%write a mesh function: this is a one element mesh now for HW6

t = 0.01;

Nel = size(e,1);%number of elements
Nnodes = size(co,1); %number of nodes
nne = 4; %number of nodes per element
dof = 2; %degree of freedom per node

%%%%%%%%%%%%PREPROCESSING END%%%%%%%%%%%
%%%Generic block: Initializes global stiffness matrix 'K' and force vector 'F'
K = sparse(Nnodes*dof,Nnodes*dof);
F = zeros(Nnodes*dof,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Assemble Global system - generic FE code
for A = 1:Nel
    coord = co(e(A,:),:);% get coord matrix for element A
    intpts = (1/sqrt(3))*[-1 -1;1 -1;1 1;-1 1]; %set of integration points, same for every element
    
    %local stiffness matrix and force vector 
    F_body = [0;0];
    localbodyf = localbody(intpts,coord,F_body,t); %hw6, p2 
    F_tr = [-11439600;0];
    tr_node = find(coord(:,1)==0.904);%find nodes (for which x = 2, HW6) where traction is applied
    localtracf = localtraction(tr_node,coord,F_tr,t); % HW6, p3
    localforce = localbodyf + localtracf; 
    
    localstiffness = localstiffnessmat(intpts,coord,t); %HW6, p4
    
    
    %DONT TOUCH BELOW BLOCK!! Assembles the global stiffness matrix, Generic 
    %block which works for any element

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
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%BOUNDARY CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%

tr_node = find(co(:,1)==0); %find nodes (for which x = 0) that need to be fixed
deletedofs = [2*tr_node-1;2*tr_node];%for these nodes, list of degrees of freedom
K(deletedofs,:) = [];
K(:,deletedofs) = [];
F(deletedofs,:) = [];
%%%%%%%%%%%%%%%%%%%BOUNDARY CONDITIONS END%%%%%%%%%%%%%%%%%%%%%%%

%solve for displacement unknowns (uk)
uk = K\F;
%expand u to include deleted displacement bcs
u = ones(Nnodes*dof,1);
u(deletedofs) = 0;
I = find(u == 1);
u(I) = uk;


figure()

hold on
patch('Faces',e,'Vertices',co,'FaceColor','none')
nodes_displaced = co + reshape(u,2,[])';
patch('Faces',e,'Vertices',nodes_displaced,'FaceColor','none','EdgeColor','Red')
hold off
axis equal
%%%Post-Process for stress

E = 70e9;
nu = 0.3;
D = E/(1-nu^2) * [ 1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2 ];

figure()
hold on
axis equal
for i = 1:Nel
    int_Coord = [-1,-1;
                 1,-1;
                 1,1;
                 -1,1];
    quad_coord = nodes_displaced(e(i,:),:);

    tau_vec = zeros(1,4);
    for j = 1:size(int_Coord,1)
        [~,~,B] = element(int_Coord(j,1),int_Coord(j,2),quad_coord);

        
        u_node = [2*e(i,:)-1;2*e(i,:)];
        q = u(u_node(:));
        tau = D*B*q;
        
        tau_von_mise = sqrt(tau(1).^2 - (tau(1).*tau(2)) + tau(2).^2 + (3.*tau(3).^2));

        tau_vec(j) = tau_von_mise;
    end
    
    patch('Faces',[1,2,3,4],'Vertices',quad_coord,'FaceVertexCData',tau_vec','FaceColor','interp','EdgeColor','none')
end
colorbar

%% Functions Declared
function [N, J, B] = element(xi, eta, coords) %hw6, p1
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

function f = localbody(intpts,coords,F_body,t)
    f = zeros(8,1);

    for i = 1:size(intpts,1)
        [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
        f = f + N'*F_body*det(J)*t;
    end
end

function k = localstiffnessmat(intpts,coords,t)
    k = zeros(8,8);
    E = 70e9;
    nu = 0.3;
    D = E/(1-nu^2) * [ 1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2 ];
    for i = 1:size(intpts,1)
        [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
        k = k + B'*D*B*det(J)*t;
    end
end

function f = localtraction(tr_node,coords,F_tr,t)
    f = zeros(8,1);
    if length(tr_node) > 0 %if nodes have traction
        intpts = [1 1/sqrt(3);1 -1/sqrt(3)]; 
        for i = 1:size(intpts,1)
            [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
            detJstar = sqrt(J(2,1)^2 + J(2,2)^2); 
            f = f + N'*F_tr*detJstar*t;
        end
    end
end



% function [co,e,elem,loc] = genmesh(nelem)
%     %for HW7, you need to write a meshing code
%     %make sure the elements are numbered ccw such that the face on which traction 
%     is applied has the second and third nodes of the element
%     %elem and loc are the element number containing (1,0) and the integration
%     %point at that location 
%     co = [ 0 0
%            2 0
%            0 2
%            2 1];
%     e = [ 1 2 4 3];
%     elem = 1;
%     loc = [0 -1];
% end