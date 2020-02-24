function [P, L, U] = PLU(A)
    % PLU Factoritza una matriu A en una triangular superior i una inferior. Retorna un vector P en que a cada posicio hi ha la fila que shi ha de posar. Servira per permutar els coeficients indepents.

    %Primer sabem la dimensio de la matriu
    [dim, n] = size(A);

    if dim ~= n
        error('PLU NECESSITA UNA MATRIU QUADRADA. Revisar els arguments');
    end

    U = A; L = eye(dim); P = [1:dim];

    %Aqui fem apareixer una columna de zeros en la columna i.

    for i = 1:1:dim

        %Abans de fer gauss, la fila amb l'exponent mes alt a dalt i recordarem les rotacions en el vector P

        [maxValue, maxIndex] = max(abs(U(i:end, i))); %maxValue es el valor maxim, i maxIndex es la posicio del valor

        maxIndex = maxIndex + i - 1; %Com que nomes mirem el maxim desde la diagonal
        %cap a baix, hem de fer aquestes sumes a maxIndex per obtenir la posicio en el vector columna
        P(i) = maxIndex;
        U([i maxIndex], :) = U([maxIndex i], :); %fila maxIndex passa a ser la fila i, i viseversa
        L([i maxIndex], :) = L([maxIndex i], :);

        for j = i + 1:1:dim%j=files
            %Per cada fila (de 2 a dim) fem apareixer un zero
            %Primer element de la fila
            L(j, i) = U(j, i) / U(i, i);
            U(j, :) = U(j, :) - (U(j, i) / U(i, i)) * U(i, :);

        end

    end

end
