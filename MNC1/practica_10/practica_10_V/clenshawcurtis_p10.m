function V = clenshawcurtis_p10(a, b, N, fx)
    %{ 
    Cuadratura de Clenshaw-Curtis aplicada al problema: No li passem
    la funció com a input sinó directament el valor de la funció
    en els nodes. Això ens permetrà no haver de reescriure la 
    funció, que depèn de x, y i z, cada cop que vulguem variar
    aquests valors.
    %}
Wj = [];
p = (1/((N^2)-1));
for j=0:N
    if j==0||j==N
        Wj = [Wj p];
    else
        suma = 0;
        % n serà  sempre par per aixo podem dividir per 2.
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