function x = qrSolve(UR, b)
    [n, m] = size(UR);

    if n > m
        R = triu(UR(1:m, :)); %seleccionem la matriu R barret
        M = m;
    else
        R = triu(UR);
        M = m-1;
    end
    
    for ii=1:M  %ens fixem en el nombre de columnes que sera m, 
        u = [1; UR(ii+1:end, ii)];
        gama = 2/(norm(u))^2;
        vAux = gama * u' * b(ii:end);
        b(ii:end) = b(ii:end) - u * vAux;
    end
    
    x= BS(R,b(1:m));
    
end


% S'ha d'introduir el vector b en format columna