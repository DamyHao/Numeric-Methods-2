clear all;
clc;
close all;

q = 2;
N = 24;
k = [-N/2: N/2 - 1];
L = diag(k.^2);

for ii = 1:1:N - 2
    L(ii+2, ii) = q;
end


for jj = 1:1:N - 2
    L(jj, jj+2) = q;
end

[VEC, VAL] = eig(L);
val = diag(VAL);

jj = [0:N-1]';
x = 2*pi*jj/N;

PHI = zeros(N);
icol = 0;
for kk = [-N/2:1:N/2-1]
    icol = icol +1;
    PHI(:, icol) = exp(1i*kk*x);
end
figure
phi1 = PHI*VEC(:,2); phi1 = phi1/max(phi1);
plot(x, phi1);

figure;
for t= 1:1:3
    plot(x, DFS(VEC(:,t),x)/max(DFS(VEC(:,t),x)));
    hold on;
end
hold off