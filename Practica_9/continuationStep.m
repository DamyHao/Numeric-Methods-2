function [y, iconv] = continuationStep(fun, y0, y1, s, tol, itmax)

%NORMA MODIFICADA PER LA PRACTICA 9!!!

    it = 1;
    tolk = 1;
    v = y1 - y0;
    yp = y1 + v * s; % Si s = 1 conseguim que la separaci� entre solucions sigui el maxim de "constant"
    xk = yp;
    XK = [];
    a=0; b=1;
    % A part de les ecuacions que teniem en el nnewton normal, li
    % imposareem que el preoducte escalar entre v i (xk(punt
    % buscat)- xk(predictor)) sigui 0
    n = length(y0)+1;
    normSpecial=@(u)(sqrt(cuadratura_cc(a, b, n-2, u.^2)));
    
    while it < itmax && tolk > tol
        J = jaco(fun, xk); % Jacobia en la posicio anterior

        J = [J; v'];
        fk = [fun(xk); v' * (xk - yp)]; % TODO: Copiat de teoria
        [P, L, U] = PLU(J);
        Dx = pluSolve(L, U, P, -fk); %Solucio de la ecuacio J*Dx = -fk
        
        %DxM = J\-fk;
        xk = xk + Dx;
        XK = [XK, xk];

        tolk = normSpecial(Dx); % Mirem la distancia entre el anterior i l'actual
        it = it + 1;
    end

    y = xk;

    %Retornem si convergeix o no per modificar la s si cal:
    if it <= itmax && tolk < tol
        iconv = 0; %OK
    else
        iconv = 1; %No em arribat a enlloc
    end

end
