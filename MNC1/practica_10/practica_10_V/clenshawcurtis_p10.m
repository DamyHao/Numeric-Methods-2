function V = clenshawcurtis_p10(a, b, N, fx)
    %{ 
    Cuadratura de Clenshaw-Curtis aplicada al problema: No li passem
    la funci� com a input sin� directament el valor de la funci�
    en els nodes. Aix� ens permetr� no haver de reescriure la 
    funci�, que dep�n de x, y i z, cada cop que vulguem variar
    aquests valors.
    %}
Wj = [];
p = (1/((N^2)-1));
for j=0:N
    if j==0||j==N
        Wj = [Wj p];
    else
        suma = 0;
        % n ser� sempre par per aixo podem dividir per 2.
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

V = Wj*fx';
V = ((b-a)/2)*V;

end