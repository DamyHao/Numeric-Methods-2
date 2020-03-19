%% Damia Casas Casajuana P1
clear all;
C = [;];
% Mida de la matriu quadrada:
% En aquest cas poso 8 perque ho diu l'enunciat
N = 8;

% Ara anire omplint la matriu per cada posicio:
for n = 0:1:N
    for k = 0:1:N
        % Aplico la definció de la Matriu per l'enunciat:
        if n == k || k == 0
            C(n+1,k+1) = 1;
            %disp('1');
        elseif k > n
            C(n+1,k+1) = 0;
        else
            C(n+1,k+1) = C(n,k) + C(n, k+1);
        end
    end
end
C