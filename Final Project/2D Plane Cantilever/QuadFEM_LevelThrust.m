clc
clear
close all
%% New 
hw7
function hw7 %function to do convergence study
        stress_max = [];
        for nelem = 2:2:50
            [max_stress] = FEMquad(nelem,50);
            stress_max = [stress_max,max_stress];
        end
        figure()
        plot((2:2:50).^2,stress_max);
        legend('\tau_x','\tau_y','\tau_{xy}')
        xlabel('Number of Elements')
        ylabel('Stress Unit (Pa)')

        figure(3)
        subplot(3,1,1)
        title('\tau_x (Pa)')
        colorbar
        subplot(3,1,2)
        title('\tau_y (Pa)')
        ylabel('y(m)')
        colorbar
        subplot(3,1,3)
        title('\tau_{xy} (Pa)')
        colorbar
        xlabel('x(m)')
end

function [max_stress] = FEMquad(nelem,ENDELEM)
    %Main FEA code for hw6, input nelem is the number of elements along x-axis
    %only even nelem here
    [co,e] = genmesh(nelem);

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
        
        tr_node = find(coord(:,1)==max(co(:,1))/2);%find nodes (for which x = 2, HW6) wheretraction is applied
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
    tr_node_1 = find(co(:,1)==0); %find nodes (for which x = 0) that need to be fixed
    tr_node_2 = find(co(:,1)==1);
    tr_node = [tr_node_1;tr_node_2];
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
    D = E/(1-nu^2).*[1, nu, 0;nu,1,0;0,0,(1-nu)/2];
    
    if nelem == ENDELEM
        figure()
        axis equal
        hold on
        patch('Faces',e,'Vertices',co,'FaceColor','none','EdgeColor','k')
        patch('Faces',e,'Vertices',co + reshape(u.*1E6,2,[])','FaceColor','none','EdgeColor','r')
        hold off
        title('Resulting Deformation (Magnified x1E6)')
        xlabel('x')
        ylabel('y(1/1E6 m)')
        xlabel('x(1/1E6 m)')
    end

    if nelem == ENDELEM
        figure()
        hold on
        colorbar
        axis equal
    end


    
    stress_vec = []; %Stress at all nodes of the bar

    for i = 1:size(e,1)
        a = e(i,:); %nodes in element containing 1,0
        dofs = [2*a-1;2*a];
        q = u(dofs(:));
        loc = [-1,-1;
                1,-1;
                1, 1;
               -1, 1];
        stress_local = zeros(3,4); %In CCW Element ordering per element
        
        for j = 1:4
            [N, J, B] = element(loc(j,1), loc(j,2), co(a,:)); %find B at (1,0)
            stress = [D*B*q]'; %find stress and store in a matrix indexed against number of elements in ele matrix
            stress_local(:,j) = stress;
        end
        
        stress_vec = [stress_vec, stress_local];

        if nelem == ENDELEM
            subplot(3,1,1)
            patch('Faces',[1,2,3,4],'Vertices',co(a,:),'FaceVertexCData',stress_local(1,:)','FaceColor','interp','EdgeColor','none')
            
            subplot(3,1,2)
            patch('Faces',[1,2,3,4],'Vertices',co(a,:),'FaceVertexCData',stress_local(2,:)','FaceColor','interp','EdgeColor','none')
         
            subplot(3,1,3)
            patch('Faces',[1,2,3,4],'Vertices',co(a,:),'FaceVertexCData',stress_local(3,:)','FaceColor','interp','EdgeColor','none')

        end
    end

    hold off
    
    [~,tau_x_max_indx] = max(stress_vec(1,:));
    [~,tau_y_max_indx] = max(stress_vec(2,:));
    [~,tau_xy_max_indx] = max(stress_vec(3,:));

    max_stress = [stress_vec(1,tau_x_max_indx);stress_vec(2,tau_y_max_indx);stress_vec(3,tau_xy_max_indx)];

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
    t = 0.01;
    a = 5;
    f_body = 0;
    for i = 1:size(intpts,1)
        [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
        f = f + N'*[0;f_body]*det(J)*t;
    end
end

function k = localstiffnessmat(intpts,coords) %hw6 p4
    k = zeros(8,8);
    t = 0.01;
    E = 70e9;nu = 0.3;
    D = E/(1-nu^2).*[1, nu, 0;nu,1,0;0,0,(1-nu)/2];
        for i = 1:size(intpts,1)
        [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
        k = k + B'*D*B*det(J)*t;
        end
end

function f = localtraction(tr_node,coords) %hw6 p3
f = zeros(8,1);
    if length(tr_node) > 0 %if nodes have traction
    
    t = 0.01;
    intpts = [1 1/sqrt(3);1 -1/sqrt(3)];
        for i = 1:size(intpts,1)
            [N, J, B] = element(intpts(i,1), intpts(i,2), coords);
            detJstar = sqrt(J(2,1)^2 + J(2,2)^2);
            f = f + N'*[0;1E3]*detJstar*t;
        end
    end
end

function [co,e] = genmesh(nelem)
    [co,e] = buildMesh(nelem);
end