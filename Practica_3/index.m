%% PR�CTICA 3: GRAM-SCHMIDT AND REORTHOGONALIZATION

%H=hilb(n); %Hilbert matrice
%k=cond(hilb(n); %numero de condicionament de H
% surfc comanda curiosa mantlab

format long g

n = 6;


for n = 2:1:12
    hilbertt = hilb(n);
    [S, Q, R] = cgs(hilbertt);
    disp(norm(eye(n) - Q' * Q));
end




