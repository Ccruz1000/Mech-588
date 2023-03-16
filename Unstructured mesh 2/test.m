clc;
clear all;
close all;
%% Reading Mesh file

meshfile_name = sprintf('Face-Cell/mech511-square-verycoarse.mesh');     % reading mesh file name

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
connection=connection;
k=0;
l=0;
for i=1:num_edge
    ind=i-1;
    if(connection(i,2)==-1)%-1)
        k=k+1;
        bdry_edge(k,:)=connection(i,:);
    else
        l=l+1;
        interior_edge(l,:)=connection(i,:);  
    end
end



for i=1:num_interior
    
    ind_edge_cell(i,1)=interior_edge(i,1);
    ind_edge_cell(i,2)=interior_edge(i,2);
    int_edg_vert(i,1)=interior_edge(i,3);
    int_edg_vert(i,2)=interior_edge(i,4);
end

for i=1:num_bcedge
    ind_edge_cell_bc(i,1)=bdry_edge(i,1);
    ind_edge_cell_bc(i,2)=bdry_edge(i,2);
    int_bc_vert(i,1)=bdry_edge(i,3);
    int_bc_vert(i,2)=bdry_edge(i,4);
end

ind_edge_cell=ind_edge_cell+1;

neigh_cells=zeros(num_cells,3);

for i=1:num_interior

    i1=ind_edge_cell(i,1);
     
    i2=ind_edge_cell(i,2);

    if (neigh_cells(i1,1)==0)
    neigh_cells(i1,1)=i2;
    elseif(neigh_cells(i1,1)~=0 & neigh_cells(i1,2)==0)
    neigh_cells(i1,2)=i2;
    else
        neigh_cells(i1,3)=i2;
    end

    if (neigh_cells(i2,1)==0)
    neigh_cells(i2,1)=i1;
    elseif(neigh_cells(i2,1)~=0 & neigh_cells(i2,2)==0)
    neigh_cells(i2,2)=i1;
    else
        neigh_cells(i2,3)=i1;
    end

end

for i=1:num_cells
    if(neigh_cells(i,1)==0)
        neigh_cells(i,1)=-1;
    end
    if(neigh_cells(i,2)==0)
        neigh_cells(i,2)=-1;
    end
    if(neigh_cells(i,3)==0)
        neigh_cells(i,3)=-1;
    end
end

for i=1:num_cells

    i1=i-1;
    z=0;
    for j=1:3
        i2=neigh_cells(i,j)-1;
        
         for k=1:num_edge
             
             if(connection(k,1)==i1 & connection(k,2)==i2)
                 z=z+1;
                 cell_edge(i,z)=k;
             elseif(connection(k,1)== i2 & connection(k,2)==i1 )
                 z=z+1;
                 cell_edge(i,z)=k;
             end
         end
    end
end

for i=1:num_cells
    Max_edge=max(cell_edge(i,:));%+1;
    Min_edge=min(cell_edge(i,:));%+1;
    
    cell_vert(i,1)=connection(Max_edge,3);
    cell_vert(i,2)=connection(Max_edge,4);
    cell_vert(i,3)=connection(Min_edge,3);

end





