%% PRACTICA 3: GRAM-SCHMIDT AND REORTHOGONALIZATION

clear all
close all
format long g

%% A) 

error = [];
N = 2:1:12;
for n = N
    hilbertt = hilb(n); %fem la matriu de Hilbertt
    [S, R, Q] = cgs(hilbertt); 
    
    error = [error, norm(eye(n) - Q' * Q)];
end

figure(1)
semilogy(N,error)
hold on;

%Producte Q'Q hauria de donar la identitat perq els vectors columna son una
%base ortonormal.
%No és aixi ja que els vectors de la matriu de hilbert es gran i densa i es
%probable que siguin practimant paralels

%% B)

error=[];

for n = N
    hilbertt = hilb(n); %fem la matriu de Hilbertt
    [Q, R] = gsr(hilbertt); 
    
    error = [error, norm(eye(n) - Q' * Q)];
end

figure(1)
semilogy(N,error)
hold on;

%% C)

error=[];

for n = N
    hilbertt = hilb(n); %fem la matriu de Hilbertt
    [Q, R] = qr (hilbertt); 
    
    error = [error, norm(eye(n) - Q' * Q)];
end

figure(1)
semilogy(N,error)
hold on;

title('REPRESENTACIO DELS ERRORS DELS DIFERENTS METDOES');
xlabel('nombre de files de la matriu');
ylabel("norma de la diferencia identitat-Q'Q");
legend('cgs','gsr','qr');





