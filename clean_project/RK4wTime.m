function V = RK4wTime(vn0, h, fun, desiredPoints)
    % Algoritme per resoldre ODEs de PVI.
    % En aquest cas tornara els calculs corresponents al temps
    %  plot(V(1,:), V(2, :));
    %  |--|--|--|--|...   "|": points  "--": steps
    %    => desiredPoints - 1 = steps (en aquest cas RK1 o AB1)
    %   Tfinal sera vn(0) + h*(desiredPoints-1)
    % Inputs:
    %   vn0: introduim vn0 columna i en la posicio 0 el instant 0
    %   h: increment de temps. Estara equiespaiat
    %   fun: Funcio f que dona la derivada: dv/dt = f(t, v(t))
    %   desiredPoints: nombre de punts que treura (comptant el que ja li
    %   donem). Es en realitat desiredPoints = steps+1 on steps son els pasos que fara.
    % Outputs:
    %   v: matriu amb els punts com a columnes EN LA PRIMERA FILA HI HAURA EL TEMPS
    % error del hordre de O(h^4) Es pot comprovar fent plot log log i
    % mirant pendent

    % Prellocating memory to gain speed
    V = zeros(length(vn0), floor(desiredPoints));
    V(:, 1) = vn0;
    %v = [vn0];
    vn = vn0(2:end);
    t = vn0(1);

    for ii = 1:desiredPoints - 1% Si li demanem un punt fara 0 iteracions
        a = h * fun(vn);
        b = h * fun(vn +a / 2);
        c = h * fun(vn + b / 2);
        d = h * fun(vn + c); % 4 avaluacions de fun

        vn1 = vn + (1/6) .* (a + 2 * b + 2 * c + d);
        V(2:end, ii + 1) = vn1; %cada columna es un "punt"
        V(1, ii + 1) = t; %Hi posem el temps
        t = t + h;
        V(1, ii + 1) = t; %Hi posem el temps

        vn = vn1;
    end

end
