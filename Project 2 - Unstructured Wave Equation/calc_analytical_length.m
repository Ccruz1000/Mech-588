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
normal(i, 1) = (1 / Edge_length(i)) * (Vertex_coord(Edge_vertex(i, 1), 2) - Vertex_coord(Edge_vertex(i, 2), 2));
normal(i, 2) = (1 / Edge_length(i)) * (-Vertex_coord(Edge_vertex(i, 1), 1) + Vertex_coord(Edge_vertex(i, 2), 1));
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