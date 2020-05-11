function vn1 = gravFunctionV(vn)
    % (r1, r2, v1, v2)
    % Funcio gravetat per la practica 10:
    m1 = 1; m2 = 2; m3 = 1;
    x = vn(1); y = vn(2);
    vtn = m1 / (x^2 + (1 - y)^2).^(3/2).* [-x, 1 - y] + m2 / (x^2 + (1 + y)^2).^(3/2) .* [-x, -1 - y] + m3 / (y^2 + (1 - x)^2)^(3/2) .* [1 - x, -y];
    vn1 = [vn(3); vn(4); vtn(1); vtn(2)];
end
