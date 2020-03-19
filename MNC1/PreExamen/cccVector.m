function integral = cccVector(a, b, funxk)
%   Clenshaw-Curtis pero entregant un vector (si es de Chevyshev millor jeje):

N = length(funxk) - 1;
Wj = [];
p = (1/((N^2)-1));
for j=0:N
    if j==0||j==N
        Wj = [Wj p];
    else
        suma = 0;
        % n serà sempre par per aixo podem dividir per 2.
        for k=0:N/2
            if k == 0 || k == N/2
                suma = suma + (1/2)*(1/(1-4*(k^2)))*cos((2*pi*k*j)/N);
            else 
                suma = suma + (1/(1-4*(k^2)))*cos((2*pi*k*j)/N);
            end
            
        end
        Wj= [Wj (4/N)*suma];
    end
end

integral = Wj*funxk';
integral = ((b-a)/2)*integral;

end

