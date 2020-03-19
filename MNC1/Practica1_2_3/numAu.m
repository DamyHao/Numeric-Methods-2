function [num_Au, iter] = numAu(tol, max_iterations)
n1 = 1;
n2 = 1;
it = 1;
result1 = 1;
result2 = 1;
aux = 0;
continuar = 0;

while it < max_iterations && continuar == 0
    aux = n1 + n2;
    n1 = n2;
    n2 = aux;
    result2 = n2/n1;
    
    %disp('Resultat 1:')
    %disp(result1);
    %disp('Resultat 2:')
    %disp(result2);
    %Comprovar si estem a tolerancia mínima. 
    if abs(result2 - result1) < tol
        continuar = 1;
    end
         
    %Actualitzar comptadors:
    it = it + 1;
    result1 = result2;
end
iter = it;
num_Au = result2;
end