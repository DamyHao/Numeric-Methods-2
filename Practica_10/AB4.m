function V = AB4(vn, h, fun, desiredPoints)
    % vn = [vn, vn+1, vn+2, vn+3] matriu de vectors columna (els 4 punts anteriors al desitjat)!!!!
    % necesitem 4 resultats anteriors:
    % ELS 4 RESULTATS HAN DETAR SEPARATS H
    % Utilitzem els coeficients ja calculats a clase:
    % b0 = -3/8; b1 = 37/24 ; b2 = -59/24; b3 = 55/24;
    betas = [-3/8, 37/24, -59/24, 55/24];
    V = vn;

    for i = 1:desiredPoints-4
        vNext = V(:, end) + h .* (betas(1) .* fun(V(:, end - 3)) + betas(2) .* fun(V(:, end - 2)) + betas(3) .* fun(V(:, end - 1)) + betas(4) .* fun(V(:, end)));
        V = [V, vNext]; %cada columna es un punt (x, y,vx,vy) i per tant V sera una matriu de 4 files
    end

end
