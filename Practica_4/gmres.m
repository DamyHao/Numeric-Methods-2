function [x, k] = gmres(Afun, V, tol, dimkyl)
    disp(V)
    normaB = norm(b); 
    res = 1;
    k = 1; 
    Q = [];
    Q(:,1) = b/norm(b); 
    disp(Q)
    H = [];
    hkp1 = 1;

    while res > tol & dimkyl > k & hkp1 > eps
        disp("finsaqui");
        disp(Q(:, k))
        Q(:, k + 1) = Afun(Q(:, k));
        
        % feval(Afun(), Q(:, k)); %primera versio del vector k+1

        %iniciem gsr

        h = Q(:, 1:k)' * Q(:, k + 1); %producte escarlar, s'ha de transposar per a que coincideixin les dimensions

        Q(:, k + 1) = Q(:, k + 1) - Q(:, 1:k) * h%segona versi� del vector k+1

        %iniciem la reorth.

        s = Q(:, 1:k)' * Q(:, k + 1); %no es igual que l'h ja que el segon terme ha canviat (agafem segona versio del vector k+1)

        Q(:, k + 1) = Q(:, k + 1) - Q(:, 1:k) * s;

        h = h + s;

        % creem la matriu de Hessemberg (m+1 x m)

        hkp1 = norm(Q(:, 1 + k));

        H(1:k + 1, k) = [h; hkp1]; % La posicio per sota la diagonal de la columna k indica la norma del vector q_k+1

        Q(:, k + 1) = Q(:, k + 1) / hkp1;

        e1b = [normaB; zeros(k, 1)]; % A la funcio zeros li donem (files, columnes). Per tant estem generant una matriu.

        UR = qrFact(H(1:k + 1, 1:k)); % no hi ha res mes pero l'alavaro ho fica

        y = qrSolve(UR, e1b);

        x = Q(:, 1:k) * y;

        k = k + 1;

        % Serveix per avaluar
        res = norm(feval(Afun, x) - b);

    end

end
