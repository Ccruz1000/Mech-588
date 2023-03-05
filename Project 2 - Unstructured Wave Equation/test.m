clc;
clear all;
close all;
%% Reading Mesh file

meshfile_name = sprintf('mech511-square-verycoarse.mesh');     % reading mesh file name

fid = fopen(meshfile_name,'r');

cur_line=fgets(fid);

results= sscanf(cur_line,'%d %d %d %d')';    % Reading the number nodes line

num_cells=results(1);

number_nodes=results(4);

num_bcedge=results(3);

num_edge=results(2);       % all

num_interior=num_edge-num_bcedge;


vertices = fscanf(fid,'%f',[2 number_nodes])';

%num_connect=num_edge+num_bcedge;

connection = fscanf(fid,'%f',[4 num_edge])';
k=0;
l=0;
for i=1:num_edge
    ind=i-1;
    if(connection(i,2,:,:)==-1)
        k=k+1;
        bdry_edge(k,:,:,:)=connection(i,:,:,:);
    else
        l=l+1;
        interior_edge(l,:,:,:)=connection(i,:,:,:);  
    end
end

cell_con=zeros(num_interior,4);

for i=1:num_interior

    index=interior_edge(i);

    cell_con(i,1)=interior_edge(i);
    cell_con(i,2)=interior_edge(i,2);

    k=2;

    for j=i+1:num_interior

        if(interior_edge(j)==index)

            interior_edge(j)=nan;

            k=k+1;

            cell_con(i,k)=interior_edge(j,2);

        end
    end
    
    
end





