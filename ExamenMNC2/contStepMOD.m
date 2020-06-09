function [y, iconv] = continuationStep(fun, y0, y1, s, tol, itmax)
% Funcio que ens dona el seguent punt de una funcio vectorial amb
% parametre. Possiblement cal utiliztarla dintre un while.
% Pot ser que "reboti". Si rebota estara a la mateixa solucio
% Input: 
%   - fun:    Funcio vectorial que el ultim component de la variable independent es el parametre
%   - y0, y1: Dos punts per trobar el seguent. Han de ser columnes. La
%             ultima component ha de ser el parametre.
% Output:
%   - y:      Seguent punt. Si s = 1 la distancia sera aproximadament la
%   distancia entre y0 i y1.
%   - iconv:  0 si hem pogut trobar el seguent punt, 1 si no. En cas de que
%   sigui 1 cal modificar la s segurament
    it = 1;
    tolk = 1;
    v = y1 - y0;
    yp = y1 + v * s; % Si s = 1 conseguim que la separaci� entre solucions sigui el maxim de "constant"
    xk = yp;
    XK = [];

    % A part de les ecuacions que teniem en el nnewton normal, li
    % imposareem que el preoducte escalar entre v i (xk(punt
    % buscat)- xk(predictor)) sigui 0
    while it < itmax && tolk > tol
        J = jaco(fun, xk); % Jacobia en la posicio anterior

        J = [J; v'];
        fk = [fun(xk); v' * (xk - yp)]; % TODO: Copiat de teoria
        %[P, L, U] = PLU(J);
        %Dx = pluSolve(L, U, P, -fk); %Solucio de la ecuacio J*Dx = -fk
        Dx = J\-fk;
        xk = xk + Dx;
        XK = [XK, xk];
        tolk = norm(Dx); % Mirem la distancia entre el anterior i l'actual
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
