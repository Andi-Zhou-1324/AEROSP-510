clc
clear
close all
%% Mesh Assembly
nodes_origin = [0,0;
                2,0;
                2,1;
                0,2];

N_total = 25;
N = sqrt(N_total);

bottomEdge_x = linspace(nodes_origin(1,1),nodes_origin(2,1),N+1);
bottomEdge_y = linspace(nodes_origin(1,2),nodes_origin(2,2),N+1);

topEdge_y = linspace(nodes_origin(4,2),nodes_origin(3,2),N+1);

%Building x coordinates
coord_x = repmat(bottomEdge_x,N+1,1);
%Building y coordinates
coord_y = zeros(N+1,N+1);

for i = 1:N+1
    temp = linspace(bottomEdge_y(i),topEdge_y(i),N+1);
    coord_y(:,i) = temp;
end

node_number = 1:(N+1)^2;

node_number = reshape(node_number,N+1,N+1)';

E2N = zeros(N.^2,4);
indx = 1;
for i = 1:N
    for j = 1:N
        temp = [node_number(i,j),node_number(i+1,j),node_number(i+1,j+1),node_number(i,j+1)];
        E2N(indx,:) = temp;
        indx = indx+1;
    end
end

nodes = [reshape(coord_x',[],1),reshape(coord_y',[],1)];

patch('Faces',E2N,'Vertices',nodes,'FaceColor','none')

%% Force Assembly

