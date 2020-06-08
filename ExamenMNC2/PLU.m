function [P, L, U] = PLU(A)
    % PLU Factoritza una matriu A en una triangular superior i una inferior amb partial pivoting
    % INPUT:
    %   A = Matriu a factoritzar
    % OUTPUT
    %   L = Matru lower en que la diagonal son 1 i per sota hi guardem les ms.
    %       Les ms son els coeficients multiplicats per la primera columna del pas
    %       enesim permeten RESTAR al primer element de la primera fila i fer una
    %       nova columna de 0s en U.
    %   U = Matriu com si li haguesim fet gauss.
    %   P = vector en que a cada posicio hi ha la fila que shi ha de posar.
    %       Servira per permutar els coeficients indepents. Si la posicio en
    %       el vector es igual al numero, significa que no hi ha hagut permutacio.
    %   Compte perque a causa de les permutacions L*U != A.

    %Primer sabem la dimensio de la matriu
    [dim, n] = size(A);

    if dim ~= n
        error('PLU NECESSITA UNA MATRIU QUADRADA. Revisar els arguments');
    end

    U = A; L = eye(dim); P = [1:dim];

    %Aqui fem apareixer una columna de zeros en la columna ii.
    for ii = 1:1:dim

        %Abans de fer gauss, la fila amb l'exponent mes alt a dalt i recordarem les rotacions en el vector P
        [maxValue, maxIndex] = max(abs(U(ii:end, ii))); %maxValue es el valor maxim, i maxIndex es la posicio del valor

        maxIndex = maxIndex + ii - 1; %Com que nomes mirem el maxim desde la diagonal
        %cap a baix, hem de fer aquestes sumes a maxIndex per obtenir la posicio en el vector columna
        P(ii) = maxIndex;
        U([ii maxIndex], :) = U([maxIndex ii], :); %fila maxIndex passa a ser la fila i, i viseversa
        L([ii maxIndex], 1:ii - 1) = L([maxIndex ii], 1:ii - 1); % En la triangular inferior no podem agafar 
        % tota la fila al complet ja que no tindra els 0 que desitgem

        for jj = ii + 1:1:dim %jj=files
            %Per cada fila (de 2 a dim) fem apareixer un zero
            %Primer element de la fila
            L(jj, ii) = U(jj, ii) / U(ii, ii); % Matru lower en que la diagonal son 1 i per sota hi guardem les ms.
            U(jj, ii:dim) = U(jj, ii:dim) - (U(jj, ii) / U(ii, ii)) * U(ii, ii:dim);

        end

    end

end