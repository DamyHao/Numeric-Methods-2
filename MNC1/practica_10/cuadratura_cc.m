function integral = cuadratura_cc(a, b, N, fun)
% Fer cuadratura Clenshaw-Curtis:

j = [0:1:N];
xcheb = cos(j.*pi./N);
xk = a + ((b-a)./2).*(xcheb+1);
fx = fun(xk);

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

integral = Wj*fx';
integral = ((b-a)/2)*integral;

end

