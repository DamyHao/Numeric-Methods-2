function x = qrSolve(UR, b)
    [n, m] = size(UR);

    if n > m
        P = UR(1:m, :);
        R = triu(P);
    else
        R = triu(UR);
    end

    vAux = gama * u' * B;
    b = b - u * vAux;
    
end
