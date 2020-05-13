function V = AB4(vn, h, fun, desiredPoints)
    % Algoritme multistep per solucionar ODEs de IVP. Més rapid que RK4.
    % Input:
    %   vn: [vn, vn+1, vn+2, vn+3] matriu de vectors columna (els 4 punts anteriors al desitjat). Calculats amb RK4 per exemple.
    %       ELS 4 RESULTATS HAN DETAR SEPARATS H.
    %   h: incrament de temps equiespaiats
    %   fun: Funcio f que dona la derivada: dv/dt = f(t, v(t)). "Part dreta de una edo de mm1"
    %   desiredPoints: nombre de punts que treura (comptant el que ja li
    %       donem).
    % Outputs:
    %   v: matriu amb els punts com a columnes
    

    % Utilitzem els coeficients ja calculats a clase:
    % b0 = -3/8; b1 = 37/24 ; b2 = -59/24; b3 = 55/24;
    betas = [-3/8; 37/24; -59/24; 55/24];
    V = vn;

    % Per evitar evaluar la fun multiples vegades farem un cua en un vector: (FIFO: first in first out)
    queue = [fun(V(:, end - 3)), fun(V(:, end - 2)), fun(V(:, end - 1)), fun(V(:, end))];

    for i = 1:desiredPoints-4
        vNext = V(:, end) + h .* (queue*betas); % Que era en realitat: (betas(1) .* fun(V(:, end - 3)) + betas(2) .* fun(V(:, end - 2)) + betas(3) .* fun(V(:, end - 1)) + betas(4) .* fun(V(:, end)))
        %Actualizem cua: el rendiment sera pobre pero es problema del matlab:
        queue(:, 1) = [];
        a = fun(vNext);
        queue(1:4, end+1) = fun(vNext);
        V = [V, vNext]; %cada columna es un punt (x, y,vx,vy) i per tant V sera una matriu de 4 files
    end
end
