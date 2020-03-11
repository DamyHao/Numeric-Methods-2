function UR = qrFact(A)
    % Factoritza la matriu A en UR, on la diagonal i per sobre es R i en els vectors u que queden per sota hi ha la informacio de les transfomracions Q aplicades.
    [n, m] = size(A);
    
    if (n > m)
        M = m; % Si no es quadrada, haurem de posar 0s sota totes les columnes
    else
        % Si factoritzem tota la matriu ens podem estalviar de posar 0s sota la ultima columna
        M = m -1;
    end

    

    for ii = 1:1:M
        B = A(ii:n, ii:m); % Es la A petita
        x = B(:, 1); % Es la primera columna de la A petita
        tau = sign(x(1)) * norm(x);
        gama = 1 + x(1) / tau; % En aquest cas gama feta com a teoria.
        u = [1; x(2:end) / (x(1) + tau)]; % Posem el punt i coma per fe que sigui columna. Posem el 1 perq el necesitem per multiplicar pero despres el treurem
        vAux = gama * u'* B;
        A(ii:n, ii:m) = B - u * vAux; A(ii+1:end,ii) = u(2:end);
    end
    
    UR = A;
