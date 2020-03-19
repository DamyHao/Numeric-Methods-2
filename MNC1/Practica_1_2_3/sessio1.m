A = [1 2 3; 1 4 2; 1 2 1];
A_t = A';
A_inv = inv(A);
A_det = det(A);
B = [A [1; 4; 8]]