clear all;
clc;
close all;
n = 26; [D, x] = chebdiffAlvaro(n); D2 = D * D; I = ones(n + 1);
c11 = 1; c12 = 0; c13 = 0; c21 = 1; c22 = 2; c23 = 0;
L = 4 * D2; M1 = -[D(1, 2:n); D(n + 1, 2:n)];
M2 = [c21 + c22 * D(1, 1), c22 * D(1, n + 1);
    c12 * D(n + 1, 1), c11 + c12 * D(n + 1, n + 1)];
M3 = [L(2:n, 1) L(2:n, n + 1)];
N = L(2:n, 2:n) + (M3 * inv(M2)) * [c22 0; 0 c12] * M1;
[EVEC, EVAL] = eig(-N); lamb = diag(EVAL);
[foo, ii] = sort(lamb); lamb = lamb(ii); EVEC = EVEC(:, ii);


p = @(x)(0 .* x + 1); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
q = @(x)(0 .* x);
r = @(x)(0 .* x);
n = 26;
a = 0; b = 1;
C = [1, 0, 0; 1, 1, 0]; %fiquem les C, p, q, r, al domini de (a,b)

[F, x, lamb1] = resoldreODEeig(n, C, p, q, r, a, b);
