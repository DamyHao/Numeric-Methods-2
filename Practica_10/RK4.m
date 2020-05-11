function v = RK4(vn0, h, fun, desiredPoints)
    % Treu la següent v. h equispaced
    % Nombre de punts que volem obtenir amb RK4
    % introduim vn0 columna
    % DesiredPoints: nombre de punts que treura (comptant el que ja li
    % donem). Es en realitat el nombre de steps-1.
    vn = vn0;
    v = [vn0];

    for i = 1:desiredPoints-1
        a = h * fun(vn);
        b = h * fun(vn +a / 2);
        c = h * fun(vn + b / 2);
        d = h * fun(vn + c);

        vn1 = vn + (1/6).* (a + 2*b + 2*c + d);
        v = [v , vn1]; %cada columna es un punt
        vn = vn1;

    end

end
