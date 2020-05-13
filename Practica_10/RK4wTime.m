function v = RK4wTime(vn0, h, fun, desiredPoints)
    % Algoritme per resoldre ODEs de PVI.
    % En aquest cas tornara els calculs corresponents al tempss
    %
    % Inputs:
    %   vn0: introduim vn0 columna SE SUPOSA QUE ES AL INSTANT 0.
    %   h: increment de temps. Estara equiespaiat
    %   fun: Funcio f que dona la derivada: dv/dt = f(t, v(t))
    %   desiredPoints: nombre de punts que treura (comptant el que ja li
    %   donem). Es en realitat desiredPoints = steps-1 on steps son els pasos que fara.
    % Outputs:
    %   v: matriu amb els punts com a columnes EN LA PRIMERA FILA HI HAURA EL TEMPS
    %
    % Truncation error: O(h³)

    % Prellocating memory to gain speed
    v = zeros(length(vn0) + 1, floor(desiredPoints));
    v(:, 1) = vn0;
    %v = [vn0];
    vn = vn0;

    t = 0;

    for i = 1:desiredPoints - 1% Si li demanem un punt fara 0 iteracions
        a = h * fun(vn);
        b = h * fun(vn +a / 2);
        c = h * fun(vn + b / 2);
        d = h * fun(vn + c); % 4 avaluacions de fun

        vn1 = vn + (1/6) .* (a + 2 * b + 2 * c + d);
        v(2:end, i + 1) = vn1; %cada columna es un "punt"
        v(1, i + 1) = t; %Hi posem el temps

        vn = vn1;

        t = t + h;
    end
end
