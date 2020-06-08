function V = ExplicitEulerwTime(vn0, h, fun, desiredPoints)
    % = RK1 o AB1 Algoritme per resoldre ODEs de PVI. 
    % V_n+1 = V_n + h*F_n
    %  plot(V(1,:), V(2, :));
    %  |--|--|--|--|...   "|": points  "--": steps
    %    => desiredPoints - 1 = steps (en aquest cas RK1 o AB1)
    %   Tfinal sera vn(0) + h*(desiredPoints-1)
    % Inputs: 
    %   t0: temps inicial (posa-hi 0 crack)
    %   vn0: introduim vn0 columna. En la primera posicio hi ha el temps.
    %   h: increment de temps. Estara equiespaiat
    %   fun: Funcio f que dona la derivada: dv/dt = f(t, v(t))
    %   desiredPoints: nombre de punts que treura (comptant el que ja li
    %   donem). Es en realitat desiredPoints = steps+1 on steps son els pasos que fara.
    % Outputs:
    %   v: matriu amb els punts com a columnes

    % Prellocating memory to gain speed
    V = zeros(length(vn0), floor(desiredPoints));
    V(:,1) = vn0;
    %v = [vn0]; 
    t = vn0(1);
    vn = vn0(2:end);
    for ii = 1:desiredPoints-1 % Si li demanem un punt fara 0 iteracions
        a = h * fun(vn);
        vn1 = vn + a;
        V(2:end, ii+1) = vn1; %cada columna es un "punt"
        t = t + h;
        V(1, ii + 1) = t; %Hi posem el temps
        vn = vn1;
    end
end