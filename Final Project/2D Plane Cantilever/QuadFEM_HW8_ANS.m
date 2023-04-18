clc
clear
close all
%% New 
hw7
function hw7 %function to do convergence study
    count = 1;
        for nelem = 2:2:50
         stress(count,:) = FEMquad(nelem);
         nel(count) = nelem;
         count = count+1;
        end
    plot(nel.*nel,stress);
end

function stress = FEMquad(nelem)
    %Main FEA code for hw6, input nelem is the number of elements along x-axis
    %only even nelem here
    [co,e,elemofinterest,loc] = genmesh(nelem);

    %write a mesh function, also
    %returns the element containing (1,0) in elemofinterest and (xi,eta) at that
    %point in loc

    Nel = size(e,1);%number of elements
    Nnodes = size(co,1); %number of nodes
    nne = 4; %number of nodes per element
    dof = 2; %degree of freedom per node

    %%%%%%%%%%%%PREPROCESSING END%%%%%%%%%%%
    %%%Generic block: Initializes global stiffness matrix 'K' and force vector 'F'
    K = zeros(Nnodes*dof,Nnodes*dof);
    F = zeros(Nnodes*dof,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Assemble Global system - generic FE code
    for A = 1:Nel
    
        coord = co(e(A,:),:);% get coord matrix for element A
        intpts = (1/sqrt(3))*[-1 -1;1 -1;1 1;-1 1]; %set of integration points, same for every element
        
        %local stiffness matrix and force vector
        localstiffness = localstiffnessmat(intpts,coord); %HW6, p3
        localforce = localbody(intpts,coord); %hw6, p2
        
        tr_node = find(coord(:,1)==2);%find nodes (for which x = 2, HW6) wheretraction is applied
        localforce = localforce + localtraction(tr_node,coord); % HW6, p4
        
        
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %%%write a code for postprocessing the stress at 1,0 and store in a
    %%%stress matrix for different number of elements
    E = 70e9;
    nu = 0.3;
    D = E/(1+nu)/(1-2*nu) * [ 1-nu nu 0 ; nu 1-nu 0 ; 0 0 1/2-nu ];
    a = e(elemofinterest,:); %nodes in element containing 1,0
    dofs = [2*a-1;2*a];
    q = u(dofs(:));
    [N, J, B] = element(loc(1), loc(2), co(a,:)); %find B at (1,0)
    stress = [D*B*q]'; %find stress and store in a matrix indexed against number of elements in ele matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
end

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
function f = localbody(intpts,coords) %hw6 p2
    f = zeros(8,1);
    t = 1;
    for i = 1:size(intpts,1)
        [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
        f = f + N'*[0;1e6]*det(J)*t;
    end
end

function k = localstiffnessmat(intpts,coords) %hw6 p4
    k = zeros(8,8);
    t = 1;
    E = 70e9;nu = 0.3;
    D = E/(1+nu)/(1-2*nu) * [ 1-nu nu 0 ; nu 1-nu 0 ; 0 0 1/2-nu ];
        for i = 1:size(intpts,1)
        [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
        k = k + B'*D*B*det(J)*t;
        end
end

function f = localtraction(tr_node,coords) %hw6 p3
f = zeros(8,1);
    if length(tr_node) > 0 %if nodes have traction
    
    t = 1;
    intpts = [1 1/sqrt(3);1 -1/sqrt(3)];
        for i = 1:size(intpts,1)
            [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
            detJstar = sqrt(J(2,1)^2 + J(2,2)^2);
            f = f + N'*[1e6;0]*detJstar*t;
        end
    end
end

function [co,e,elem,loc] = genmesh(nelem)
%generate mesh for HW7
%make sure the elements are numbered ccw such that the face on which traction is applied has the second and third nodes of the element
    [X,Y] = meshgrid(0:2/nelem:2,0:2/nelem:2);
    X = X'; Y = Y';
    co = [X(:) Y(:)];
    count = 1;
    for i = 1:(nelem+1)*nelem
        if rem(i,nelem+1) ~= 0
            e(count,:) = [i i+1 i+nelem+2 i+nelem+1];
            count = count + 1;
        end
    end
    co(:,2) = co(:,2).*(1 - co(:,1)./4);
    elem = 1+ (nelem/2);%even case, right of center
    loc = [-1,-1];%std node 1
end