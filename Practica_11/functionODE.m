function output = funtionODE(entrada)
    % Funcio de la practica 11 per resoldre la ode
    %   Input: x;y columna
    %   Output: dx/dt; dy/dt fila

    alfa = 1;
    betaC = 0.9;
    output = [entrada(2); alfa * entrada(1) + betaC * entrada(2) - entrada(1)^3 - (entrada(1)^2) * entrada(2)];

end
