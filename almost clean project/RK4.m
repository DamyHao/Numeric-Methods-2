function v = RK4(vn0, h, fun, desiredPoints)
    % Algoritme per resoldre ODEs de PVI.
    % 
    %  |--|--|--|--|...   "|": points  "--": steps
    %  t/h es steps
    % 
    % Inputs: 
    %   vn0:introduim vn0 columna
    %   h: increment de temps. Estara equiespaiat
    %   fun: Funcio f que dona la derivada: dv/dt = f(t, v(t))
    %   desiredPoints: nombre de punts que treura (comptant el que ja li
    %   donem). Es en realitat desiredPoints = steps+1 on steps son els pasos que fara.
    % Outputs:
    %   v: matriu amb els punts com a columnes

  

    % Prellocating memory to gain speed
    v = [];
    v(:,1) = vn0;
    %v = [vn0]; 
    vn = vn0;

    for i = 1:desiredPoints-2 % Si li demanem un punt fara 0 iteracions
        a = h * fun(vn);
        b = h * fun(vn +a / 2);
        c = h * fun(vn + b / 2);
        d = h * fun(vn + c);

        vn1 = vn + (1/6).* (a + 2*b + 2*c + d);
        v(:, i+1) = vn1; %cada columna es un "punt"
        vn = vn1;
    end
end
