%% P8: CASAS Y JIMENEZ
%% 1) Diferencia finita local:
%{ 
Durante las horas de pr�ctica observamos que era posible
realizar este apartado de dos maneras distintas. En primer 
lugar, aproximando la derivada con la expresi�n que aparece
en el enunciado y, en segundo lugar, mediante la construcci�n
de la matriz de derivaci�n tal y como se hizo en la pr�ctica 7.

Hemos resuelto este apartado de los dos modos.
%}
%% A partir de la f�rmula de derivaci�n aproximada:
close all;
clear all;
format long g;

comptador = 1;

figure();
for n =10:10:40
    x = linspace(-1,3,n+1);
    h = abs(x(2)-x(1));
    fd_aprox = ((12*h)^(-1))*((-1)*(fun(x+2*h))+8*fun(x+h)-8*fun(x-h)+fun(x-2*h));
    fd_aprox = fd_aprox(3:end-2);
    
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    subplot(2, 2, comptador);
    plot(x,fd);
    hold on;
    x = x(3:end-2);
    plot(x,fd_aprox);
    hold off;
    comptador = comptador + 1;
end

%% A partir de la matriz de derivaci�n:
clear all;
format rat;

comptador=1;
figure;
for n =10:10:40
    x = linspace(-1,3,n+1);
    h = abs(x(2)-x(1));
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    
    fx = fun(x);
    fx = fx';
    u = zeros(1,n+1);
    v = zeros(1,n+1);
    % No cal que toquem el primer element perque ja es 0 en els dos.
    u(2) = 8;
    u(3) = -1;
    v = -u;
    
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    subplot(2,2,comptador)
    plot(x,fd);
    hold on;
    M = toeplitz(v, u);
    fd_aprox = 1/(12*h).*M*fx;
    fd_aprox = fd_aprox(3:end-2);
    x = x(3:end-2);
    plot(x,fd_aprox);
    hold off;
    comptador = comptador+1;
end

%% Error m�ximo y representaci�n:
%{
Para resolver este apartado, repetiremos el procedimiento anterior
(solo con la matriz de derivaci�n) para una n creciente n=10:10:1000
y consideraremos, para cada n, el error m�ximo, lo acumularemos
en un vector y lo representaremos en escala doble logar�tmica en
funci�n de la n que le corresponda.
%}

clear all;
format long g;
errorsMax = [];
Ns = [10:10:1000];
figure;
for n = Ns
    x = linspace(-1,3,n+1);
    h = abs(x(2)-x(1));
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    
    %Obtenim la matriu i el vector de les derivades aproximades:
    fx = fun(x);
    fx = fx';
    u = zeros(1,n+1);
    v = zeros(1,n+1);
    u(2) = 8;
    u(3) = -1;
    v = -u;
    M = toeplitz(v, u);
    fd_aprox = 1/(12*h).*M*fx;
    % Vector de derivades exactes:
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    
    fd_aprox = fd_aprox(3:end-2);
    fd = fd(3:end-2);
    
    error = max(abs(fd_aprox'-fd));
    errorsMax = [errorsMax error];
    
end
%{
Una vez representado, y para comprovar que decae con potencia p=-4,
utilizaremos el comando regression para construirnos una recta
lo m�s semejante posible a la representaci�n del error y 
obtendremos de esa recta su pendiente.
%}
loglog(Ns, errorsMax);
title('Error m�ximo de la aproximaci�n (D. local)');
xlabel('N�mero de nodos (n)');
ylabel('Error local (E_n)');
[r,m,~] = regression(log10(Ns), log10(errorsMax))

%{
Como podemos ver, la recta generada es de buena calidad (|r|~1) y su
pendiente es m~-4, por lo que confirmamos la hip�tesis inicial: el
error decae con potencia p=-4.
%}
%% 2) Derivaci�n global:
%{
En este apartado debemos aproximar la derivada de la funci�n de
manera global, es decir, con un �nico polinomio interpolador
que contenga la informaci�n de todos los nodos (el valor de la 
funci�n en ellos) y a los que habr� que aplicar la matrix
de diferenciaci�n para un n creciente n=4:4:32.

Como hemos hecho antes, acumularemos en un vector el valor del 
error m�ximo y lo representaremos en funci�n de la n que le
corresponda.
%}
clear all;
errorsMax = [];

Ns = [4:4:32];

for n = Ns
    x = linspace(-1,3,n+1); % n+1 punts equiespaciats entre -1 i 3.
    fx = fun(x);
    fx = fx';
    
    D = diferenciacion(x);
    fd_aprox = D*fx;
    
    % Vector de derivades exactes:
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    fd_aprox = fd_aprox';
    
    error = max(abs(fd_aprox-fd));
    errorsMax = [errorsMax error];
end
figure;
loglog(Ns, errorsMax);
title('Error m�ximo de la aproximaci�n (D. global)');
xlabel('N�mero de nodos (n)');
ylabel('Error local m�ximo (E_n)');

%{
En aquesta figura �s dificil veure Runge, tanmateix el podem 
observar en el rebot de l'error que es produeix a partir de n = 28. 
Si fem el mateix proc�s per� arribant a n m�s altes (per exemple 128) 
podem veure el fenomen de Runge m�s clarament:
%}

clear all;
errorsMax = [];

Ns = [4:4:128];

for n = Ns
    x = linspace(-1,3,n+1);
    % n+1 punts equiespaciats entre -1 i 3.
    fx = fun(x);
    fx = fx';
    
    D = diferenciacion(x);
    fd_aprox = D*fx;
    
    % Vector de derivades exactes:
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    fd_aprox = fd_aprox';
    
    error = max(abs(fd_aprox-fd));
    errorsMax = [errorsMax error];
end

figure;
loglog(Ns, errorsMax);
title('Derivaci�n global para n grandes')
xlabel('N�mero de nodos (n)');
ylabel('Error local m�ximo(E_n)');

%{
En aquesta gr�fica �s f�cil veure el fenomen de Runge ja que 
apareixen errors de fins a 10^25.
%}
