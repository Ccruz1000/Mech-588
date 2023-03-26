Vertex_coord = [[3.25, 1.25]; [3, 0.25]; [2, 1]; [4.5, 0.5]; [1.5, 0]; [2.5, 2.25]];

Edge_vertex = [[6, 3]; [3, 5]; [1, 6]; [2, 4]; [3, 1]; [4, 1]; [2, 3]; [5, 2]; [2, 1]];

for i = 1:length(Edge_vertex)
x = Vertex_coord(Edge_vertex(i, 1), 1) - Vertex_coord(Edge_vertex(i, 2), 1);
y = Vertex_coord(Edge_vertex(i, 1), 2) - Vertex_coord(Edge_vertex(i, 2), 2);
Edge_length(i) = sqrt(x^2 + y^2);
end

Cell_vertex = [[1, 2, 3]; [1, 2, 4]; [1, 3, 6]; [2, 3, 5]];

for i = 1:length(Cell_vertex)
x1 = Vertex_coord(Cell_vertex(i, 1), 1);
x2 = Vertex_coord(Cell_vertex(i, 2), 1);
x3 = Vertex_coord(Cell_vertex(i, 3), 1);
Cell_centroid(i, 1) = (x1 + x2 + x3)/3;
y1 = Vertex_coord(Cell_vertex(i, 1), 2);
y2 = Vertex_coord(Cell_vertex(i, 2), 2);
y3 = Vertex_coord(Cell_vertex(i, 3), 2);
Cell_centroid(i, 2) = (y1 + y2 + y3)/3;
end

for i = 1:length(Edge_vertex)
x1 = Vertex_coord(Edge_vertex(i, 1), 1);
x2 = Vertex_coord(Edge_vertex(i, 2), 1);
Edge_centroid(i, 1) = (x1 + x2)/2;
y1 = Vertex_coord(Edge_vertex(i, 1), 2);
y2 = Vertex_coord(Edge_vertex (i, 2), 2);
Edge_centroid(i, 2) = (y1 + y2)/2;
end

for i = 1:length(Edge_length)
normal(i, 1) = (1 / Edge_length(i)) * (Vertex_coord(Edge_vertex(i, 2), 2) - Vertex_coord(Edge_vertex(i, 1), 2));
normal(i, 2) = (1 / Edge_length(i)) * (-Vertex_coord(Edge_vertex(i, 2), 1) + Vertex_coord(Edge_vertex(i, 1), 1));
end

for i = 1:length(Edge_centroid)
vel(i, 1) = Edge_centroid(i, 2) * pi;
vel(i, 2) = -Edge_centroid(i, 1) *pi;
end

for i = 1:length(Cell_vertex)
Xb = Vertex_coord(Cell_vertex(i, 2), 1);
Xa = Vertex_coord(Cell_vertex(i, 1), 1);
Yc = Vertex_coord(Cell_vertex(i, 3), 2);
Ya = Vertex_coord(Cell_vertex(i, 1), 2);
Yb = Vertex_coord(Cell_vertex(i, 2), 2);
Xc = Vertex_coord(Cell_vertex(i, 3), 1);
Cell_Area(i) = 0.5*((Xb - Xa)*(Yc - Ya) - (Yb - Ya)*(Xc-Xa));
end

for i = 1:length(Edge_centroid)
Edge_vel(i, 1) = pi*Edge_centroid(i, 2);
Edge_vel(i, 2) = -pi*Edge_centroid(i, 1);
end

Temp_cent = [100; 102; 101; 97];

% Calculate solution at edge AB (index 9), based on cell i (0) being
% upwind cell
% Difference in cell i and cell a centroids
dxa = Cell_centroid(2, 1) - Cell_centroid(1, 1);
dya = Cell_centroid(2, 2) - Cell_centroid(1, 2);

% Difference in cell i and cell b centroids
dxb = Cell_centroid(3, 1) - Cell_centroid(1, 1);
dyb = Cell_centroid(3, 2) - Cell_centroid(1, 2);

% Difference in cell i and cell c centroids
dxc = Cell_centroid(4, 1) - Cell_centroid(1, 1);
dyc = Cell_centroid(4, 2) - Cell_centroid(1, 2);

% Difference in cell averaged temperatures
dta = Temp_cent(2) - Temp_cent(1);
dtb = Temp_cent(3) - Temp_cent(1);
dtc = Temp_cent(4) - Temp_cent(1);

% Calculate values for gradient in cell i
DX2 = dxa^2 + dxb^2 + dxc^2;
DY2 = dya^2 + dyb^2 + dyc^2;
DXDY = dxa*dya + dxb * dyb + dxc * dyc;
DXDT = dxa * dta + dxb * dtb + dxc * dtc;
DYDT = dya * dta + dyb * dtb + dyc * dtc;

% Calculate gradient in cell i
gradix = (1 / (DX2 * DY2 - DXDY^2)) * (DY2 * DXDT - DXDY * DYDT);
gradiy = (1 / (DX2 * DY2 - DXDY^2)) * (-DXDY * DXDT + DX2 * DYDT);

% Calcualte solution at edge AB (index 9)
T_AB = Temp_cent(1) + gradix * (Edge_centroid(9, 1) - Cell_centroid(1, 1))+ gradiy * (Edge_centroid(9, 2) - Cell_centroid(1, 2));

% Create test for flux at edge AB
vel_test1 = [2; 0];
dot_product1 = normal(9, 1) * vel_test1(1) + normal(9, 2) * vel_test1(2);
flux_horizontal = dot_product1 * T_AB * Edge_length(9);

% Test again with velocity downward with magnitude 2\
vel_test2 = [0; -2];
dot_product2 = normal(9, 1) * vel_test2(1) + normal(9, 2) * vel_test2(2);
flux_vertical = dot_product2 * T_AB * Edge_length(9);