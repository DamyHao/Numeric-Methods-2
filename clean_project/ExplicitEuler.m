function v = ExplicitEuler(vn0, h, fun, desiredPoints)
    % = RK1 o AB1 Algoritme per resoldre ODEs de PVI. 
    % V_n+1 = V_n + h*F_n
    % 
    %  |--|--|--|--|...   "|": points  "--": steps
    % => desiredPoints - 1 = steps (en aquest cas RK1 o AB1)
    % Inputs: 
    %   vn0: introduim vn0 columna
    %   h: increment de temps. Estara equiespaiat
    %   fun: Funcio f que dona la derivada: dv/dt = f(t, v(t))
    %   desiredPoints: nombre de punts que treura (comptant el que ja li
    %   donem). Es en realitat desiredPoints = steps+1 on steps son els pasos que fara.
    % Outputs:
    %   v: matriu amb els punts com a columnes

    % Prellocating memory to gain speed
    v = zeros(length(vn0), floor(desiredPoints));
    v(:,1) = vn0;
    %v = [vn0]; 
    vn = vn0;
    for ii = 1:desiredPoints-1 % Si li demanem un punt fara 0 iteracions
        a = h * fun(vn);
        vn1 = vn + a;
        v(:, ii+1) = vn1; %cada columna es un "punt"
        vn = vn1;
    end
end