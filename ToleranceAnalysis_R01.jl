using Measurements, LinearAlgebra
° = π/180;  # Degree symbol = conversion from degree to radians

# Create a quick function to rotate vectors
vec_rotate(θ, vector) = [cos(θ) -sin(θ); sin(θ) cos(θ)] * vector

A = (30 ± 0.2)°;
B = [120, 0];
C = vec_rotate(A, [120, 0]);
D = [25, 0];
E = vec_rotate(A, [-25, 0]);
F = [0, 44.95 ± 0.05];
G = vec_rotate(A, [0, 44.95 ± 0.05]);

path1 = B + C + G + E;
path2 = F + D;

Z = path1 - path2;

norm(Z) # The length of Z