%AE 510 Class
%Author: Your instructor

clear
close all
clc

%%%%%%%%%%%%PREPROCESSING%%%%%%%%%%%
%coordinate matrix [x,y] for each node
length_bar = 0.271; %meter
height_bar = 0.01; %meter

A = length_bar.*height_bar;
alpha = 5E-6;
D_T = 100;
co = [0,0,0;
      -1.47,0,0;
      0.34,-0.71,-0.55;
      -0.73,-0.71,-0.55;
      -1.81,-0.71,-0.55;
      -0.24,-1.42,-1.09;
      -1.22,-1.42,-1.09;];

figure()
scatter3(co(:,1),co(:,2),co(:,3));

E2N = [1,2;
       1,3;
       1,4;
       2,5;
       2,4;
       4,5;
       3,4;
       5,7;
       7,6;
       3,6;
       4,7;
       4,6];

figure
plot_edge(co,E2N,'k',true);
view(3)
E = 70E9;

%element-node connectivity matrix (and area for each truss in column 3)
e = [E2N,repmat(A,size(E2N,1),1)];

Nel = size(e,1);%number of elements
Nnodes = size(co,1); %number of nodes
nne = 2; %number of nodes per element
dof = size(co,2); %degree of freedom per node

%%%%%%%%%%%%PREPROCESSING END%%%%%%%%%%%

%%%Generic block: Initializes global stiffness matrix 'K' and force vector 'F'
K = zeros(Nnodes*dof,Nnodes*dof);
F = zeros(Nnodes*dof,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%Assemble Global system - generic FE code for 2D and 3D trusses
for A = 1:Nel
 
    n = (co(e(A,2),:) - co(e(A,1),:));   %n = [x2-x1;y2-y1] 

    L = norm(n); %length of truss 
    
    n = n./L; %n = [sin,cos]    
    Area = e(A,3); %area
    
    k11 = (E*Area/L)*(n'*n);%k matrix part = EA/L*[c^2 cs;sc s^2]
     
    %local stiffness matrix and force vector
    localstiffness = [k11 -k11;-k11 k11];    %full local stiffness matrix
    localforce = zeros(nne*dof,1);%external forces are added at the end, so leave as zeros. If temp changes, modify for thermal expansion 
    %     localforce = localforce + ((E.*Area*alpha*D_T).*[-n./L n./L])';
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
    
end

K_copy = K;
%%%%%%%%%%%%%%%%%%%BOUNDARY CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%

Weight_Drone = 43.0913;
Thrust = 3.*Weight_Drone.*9.81;
Thrust_per_motor = Thrust./6;

%external forces
force_node = [4];


F = assignForce(F,force_node, repmat([0,0,-Weight_Drone*9.81],size(force_node,1),1),dof);

%Apply displacement BC by eliminating rows and columns of nodes 3-4 (corresponding to
%degrees of freedom 5 to 8) - alternative (and more generic method) is the penalty approach, or
%static condensation approach - see later class notes

Nodes_fixed = [1;2;3;5;6;7];
deletedofs = fixNode(Nodes_fixed,dof);

K(deletedofs,:) = [];
K(:,deletedofs) = [];
F(deletedofs,:) = [];

%%%%%%%%%%%%%%%%%%%BOUNDARY CONDITIONS END%%%%%%%%%%%%%%%%%%%%%%%

 
%solve for displacement unknowns (uk)
uk = K\F

%expand u to include deleted displacement bcs
u = ones(Nnodes*dof,1);
u(deletedofs) = 0;
I = find(u == 1);
u(I) = uk;

u_copy = [uk;zeros(4,1)];
%%%%%%%%%%%%%POSTPROCESSING%%%%%%%%%%%%%%%%%%%%%

%%%Step 6:Postprocess results
for i = 1:Nel

    %get data about truss i
    n = (co(e(i,2),:) - co(e(i,1),:));    
    L = norm(n);    
    n = n./L;
    
    Area = e(i,3);
    n1 = e(i,1);n2 = e(i,2);%global numbers for node 1 and 2 of truss i

    d = [u(n1*dof-dof+1:n1*dof)' u(n2*dof-dof+1:n2*dof)']';%displacements of the two nodes
    sigma(i) = E*([-n./L n./L]*d) - E*Area*alpha*D_T;%stress formula, If temp changes, modify for thermal expansion 
end
sigma
hold off
%%% Plotting Results
figure
hold on
plot_edge(co,E2N,'k',true);
view(3)
plot_edge(co + reshape(u,dof,[])',E2N,'r',false);

%% Functions Declared
function [deletedofs] = fixNode(nodes,dof)
    deletedofs = [];
    for i = 1:size(nodes,1)
        n1 = nodes(i);
        temp = n1*dof-dof+1:n1*dof;
        deletedofs = [deletedofs;temp'];
    end
   
end

function [F] = assignForce(F,node,forces,dof)
    %node: an array of node number
    %forces: [nNode x dof] force array
    for i = 1:size(node,1)
        n1 = node(i);
        F(n1*dof-dof+1:n1*dof) = forces(i,:)';
    end


end

function plot_edge(nodal_list, connectivity_matrix,color,labelSwitch)
    x = nodal_list(:, 1);
    y = nodal_list(:, 2);
    z = nodal_list(:, 3);
    
    % Plot nodes with labels
    if labelSwitch
        for i = 1:size(nodal_list, 1)
            text(x(i), y(i), z(i), num2str(i), 'Color', 'b', 'FontSize', 12);
            hold on;
        end
    end

    % Plot edges
    for i = 1:size(connectivity_matrix, 1)
        start_node = nodal_list(connectivity_matrix(i, 1), :);
        end_node = nodal_list(connectivity_matrix(i, 2), :);
        plot3([start_node(1), end_node(1)], [start_node(2), end_node(2)], [start_node(3), end_node(3)], 'Color', color);
    end

    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
end
