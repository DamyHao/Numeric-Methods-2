%Exercici 13
[numAu, iteracions] = numAu(0.01, 1000);
iteracionsNecessaries = iteracions;

%Exercici 14
f = @quadrat;
g = @cub;
Resposta1 = derivadaEnPunt(f, 1);
Resposta2 = derivadaEnPunt(g, 1);