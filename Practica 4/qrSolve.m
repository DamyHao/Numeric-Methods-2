function x = qrSolve(UR, b)
    % Soluciona un sistema factoritzat amb QR. Compatible amb probleames de least squares.
    % Cal entregar-li:
    % UR (Matriu en que hi ha la R tambe les transformacions Q en les columnes decreixents on hi hauria d'haver-hi 0s)
    % b: vector dels termes independents en format columna
    [n, m] = size(UR);

    if n > m
        M = m; % La matriu no es quadrada. Servira per els least squares problem
    else
        M = m - 1; %La matriu es cuadrada i ens podem estalviar la generacio de 0s sota l'ultim element
    end

    for ii = 1:M %ens fixem en el nombre de columnes que sera m,
        u = [1; UR(ii + 1:n, ii)]; % Extreiem el vector u que conte la informacio sobre la transformacio Q
        % Apliquem les transformacions del vector u sobre el vector b
        gama = 2 / (norm(u))^2; % (Tmbe es pot generar aixi gama)
        %disp(1 + x(1) / (sign(x(1)) * norm(x)))
        vAux = gama * u' * b(ii:end);
        b(ii:n) = b(ii:n) - u * vAux;
    end

    % Finalment resolem el sistema
    R = triu(UR(1:m, 1:m)); %seleccionem la matriu R barret
    x = BS(R, b(1:m));

end
