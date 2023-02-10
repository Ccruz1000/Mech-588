# Close everything and clear
# Quadratic case to check indices
clear;
clc;
close all;

# Define function to plot the meshmetric grid being used
function plot = plot_meshmetric(neta, nxi, maxlim, minlim)
  
# Calculate y in terms of xi and eta
function y_calc = y_calc(eta, xi, neta, nxi)
y_calc = eta;
endfunction

# Calculate x in terms of xi and eta
function x_calc = x_calc(eta, xi, neta, nxi)
x_calc = xi;
endfunction

# Create vectors for xi and eta
eta = linspace(minlim, maxlim, neta);
xi = linspace(minlim, maxlim, nxi);
#eta = 0:neta; xi = 0:nxi;

# Define mesh grid for xi and eta
[ETA, XI] = meshgrid(eta, xi);

# Loop through mesh to calculate X and Y
for i = 1:(neta)
  for j = 1:(nxi)
    x(i, j) = x_calc(ETA(i, j), XI(i, j), neta, nxi);
    y(i, j) = y_calc(ETA(i, j), XI(i, j), neta, nxi);
  endfor
endfor

hold on
# Plot lines of eta
for i = 1:neta
plot(x(:,i), y(:,i))
endfor

# Plot lines of xi
for j = 1:nxi
plot(x(j,:), y(j,:))
endfor

# Write values to file
fileName = sprintf('Flux_testmesh-%dx%d.vel', neta, nxi);
fileID = fopen(fileName, 'w');
fprintf(fileID, "%d %d\n", neta, nxi);
for i = 1:(neta)
  for j = 1:(nxi)
    fprintf(fileID, "%d %d %f %f %f %f\n", i - 1, j - 1, x(i, j), y(i, j), ETA(i, j), XI(i, j));
   end
  end 
fclose(fileID);

# Store eta and xi coordinates in text file_in_loadpath
# etafileName = sprintf('ETA-XI Coordinates %dx%d.txt', neta, nxi);
#etafileId = fopen(etafileName, 'w');
#fprintf(etafileId, '%dx%d mesh', neta, nxi);
#for i = 1:neta
#  for j = 1:nxi
#    fprintf(etafileId, '%d %d %f %f\n', i-1, j-1, ETA(i, j), XI(i, j));
#   end
#  end 
#fclose(etafileId)
#printf('ETA-XI file created\n');
endfunction
