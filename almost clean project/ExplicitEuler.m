function v = ExplicitEuler(vn0, h, fun, desiredPoints)
    % Algoritme per resoldre ODEs de PVI.
    % 
    %  |--|--|--|--|...   "|": points  "--": steps
    % 
    % Inputs: 
    %   vn0:introduim vn0 columna
    %   h: increment de temps. Estara equiespaiat
    %   fun: Funcio f que dona la derivada: dv/dt = f(t, v(t)), si al
    %   funcio es vectoria ha de ser en columna
    %   desiredPoints: nombre de punts que treura (comptant el que ja li
    %   donem). Es en realitat desiredPoints = steps+1 on steps son els pasos que fara.
    % Outputs:
    %   v: matriu amb els punts com a columnes

  

    % Prellocating memory to gain speed
    v = zeros(length(vn0), floor(desiredPoints));
    v(:,1) = vn0;
    %v = [vn0]; 
    vn = vn0;
    for i = 1:desiredPoints-1 % Si li demanem un punt fara 0 iteracions
        a = h * fun(vn);
        vn1 = vn + a;
        v(:, i+1) = vn1; %cada columna es un "punt"
        vn = vn1;
    end
end
