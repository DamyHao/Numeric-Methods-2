%% P11: CASAS  Y JIMENEZ
%% 1) Integraci�n de K_theta para theta = pi/6, pi/4, pi/3, pi/3 i 3*pi/4:
clear all;
format long g;
close all;
clc;

%{
Para resolver este apartado utilizaremos tres cuadraturas distintas, 
lo cual nos permitir� dar mayor validez a nuestro resultado y
contrastar y comparar su funcionamiento para integrar:

K = @(phi)(1./((1-(k.^2).*(sin(phi)))).^(1/2));

Se resuelve la integral mediante:

- Trapezoidal Compuesta: Dado que funciona especialmente bien para
funciones peri�dicas como �sta (aparece un seno).

- Cuadratura de Clenshaw - Curtis: Dado que se fundamenta en los
polinomios de Chebyshev, m�s estables num�ricamente, y que permite as�
obtener un buen resultado para la mayor�a de casos.

- Cuadratura de F�jer: Dado que se corresponde con la cuadratura de
Clenshaw - Curtis abierta (sin incluir los extremos), y que aplicada sobre 
integrales que presentan alguna singularidad es robusta pero lenta y 
nos puede proporcionar una primera aproximaci�n al valor de la 
integral que nos sirva de referencia a la hora de aplicar otros m�todos.

Las funciones que calculan el valor de la integral mediante los
m�todos anteriormente presentados son, respectivamente, las siguientes:

function I_trap_c = trapez_comp(a, b, m, fun)
k = 0:m;
H = (b-a)/m;
xk = a+k*H;
sumatori = 0;
for i=0:m-1
    sumatori = sumatori + (fun(xk(i+1))+fun(xk(i+2)));
end
I_trap_c = (H/2)*sumatori;
end


function integral = cuadratura_cc(a, b, N, fun)
%   Fer cuadratura Clenshaw-Curtis:
%   La funcio entregada fun ha de ser capaç de tractar vectors si se li
%   donen com arguments.
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


function I_fej= fejer(a, b, N, fun)
%   Cuadratura de Fejer, que es com la de Clenshaw pero oberta.
%   La funcio entregada fun ha de ser capaç de tractar vectors si se li
%   donen com arguments.

if(1 == mod(N, 2))
    error('"N" HA DE SER PARELL!!');
end

k = 1:N;
Wj = zeros(1,N);
thk = ((2.*k-1).*pi)/(2.*N);
posicio = 0;

for j = thk
    posicio = posicio + 1;
    sumatori=0;
    for i = 1:N/2
        sumatori = sumatori + (cos(2*i*j)/(4*(i^2)-1));
    end
    
    Wj(posicio) = (2/N)*(1-2*sumatori);
end

xk = cos(thk);
phi = b + (b-a).*(1/2).*(xk-1);
f_phis = fun(phi);
I_fej = ((b-a)/2).*(Wj*f_phis');
end
%}

%L�mites de integraci�n:
a = 0;
b = pi/2;

%Requerimiento del problema: 10 cifras significativas
tol = 10^-10;

valorsTheta = [pi/6, pi/4, pi/3, pi/3, 3*pi/4];
for theta = valorsTheta
    %La funci� a integrar s'actualitza segons el valor de theta:
    k = sin(theta/2);
    K = @(phi)(1./((1-(k.^2).*(sin(phi)).^2)).^(1/2));
    
    m1 = 2;
    error = 1;
    actual = 1;
    while error > tol %TRAPEZOIDAL COMPUESTA
        I_trap_c = trapez_comp(a, b, m1, K);
        error = abs(actual - I_trap_c);
        actual = I_trap_c;
        m1 = m1+10;
    end
    
    m2 = 2;
    error = 1;
    actual = 1;
    while error > tol %CLENSHAW-CURTIS
        I_cc = cuadratura_cc(a, b, m2, K);
        error = abs(actual - I_cc);
        actual = I_cc;
        m2 = m2+10;
    end
   
    m3 = 2;
    error = 1;
    actual = 1;
    while error > tol %F�JER
        I_fej = fejer(a, b, m3, K);
        error = abs(actual - I_fej);
        actual = I_fej;
        m3 = m3+10;
    end
%Y proporciono los valores que voy obteniendo para cada valor de theta:
display = sprintf('Para theta = %d, el valor de K, seg�n el m�todo utilizado, es:', theta);
disp(display);
I_trap_c
I_cc
I_fej
end
%% 2) Aproximaci�n del periodo para peque�as oscilaciones. �ngulo a partir del cual el error de la aproximaci�n supera 5%:

% En este archivo solo se har� el proceso de encontrar el 0 mediante
% Newton. La definici�n de la funci�n a la que buscaremos el 0 es:
%{
function error = fa_funcio(theta)
% Aproximaci�n del periodo del p�ndulo para peque�as oscilaciones:
T_posc = 2*pi*((1/9.81)^(1/2));

k = sin(theta/2);
K = @(phi)(1./((1-(k.^2).*(sin(phi)).^2)).^(1/2));
a = 0;
b = pi/2;
% m son els nodes, posarem uns 200 per probar.
m = 100;
T_real = 4*((1/9.81)^(1/2))*fejer(a, b, m, K);

% Finalmente hacemos el error relativo y le restamos 0.05 porque el
enunciado pide encontrar donde se supera el 5%. As� ya preparamos la
funci�n para el metodo de Newton:
error = ((abs(T_real - T_posc)/T_real) - 0.05);
end
%}

g = 9.81;
l = 1;

k = sin(theta/2);
a = 0;
b = 2*pi;

% Visualizacion de la funcion a la que buscaremos la ra�z mediante Newton:
ys = [];
titu = [0:0.01:pi];
for i = titu
    ys = [ys fa_funcio(i)];
end
figure;
plot(titu, ys);
title('Visualizaci�n de la funci�n (error relativo(/Theta) - 0.05)');


% Hacemos la exploraci�n de Newton entre 0 y pi, sin incluir este porque
% son los valores que puede coger el angulo theta.
X = newton(0, pi-10^-5, 10^-3, 1000, @fa_funcio);
X


%% 3) Representaci�n de K_theta en funci�n del �ngulo:
%{
En este apartado deberemos calcular y representar K_theta desde theta=0
hasta valores muy proximos al �ngulo pi, para el cual prevemos que
la integral presentar� una singularidad de tipo asint�tico.

Observaremos tambi�n que el m�todo de F�jer no funciona bien en
acercarnos a la singularidad.
%}
valorsTheta = [0:0.99*pi/100:pi]; %sin llegar a valer pi.
m = 200;
K_theta_cc = zeros(1,length(valorsTheta));
K_theta_fej = zeros(1,length(valorsTheta));
n = 1;
for theta = valorsTheta
    k = sin(theta/2);
    K = @(phi)(1./((1-(k.^2).*(sin(phi)).^2)).^(1/2));
    K_theta_cc(n) = cuadratura_cc(a, b, m, K);
    K_theta_fej(n) = fejer(a, b, m, K);
    n = n + 1;
end
figure;
plot(valorsTheta, K_theta_fej, 'r');
hold on;
plot(valorsTheta, K_theta_cc, 'b');
title('Valor de la integral seg�n el �ngulo');
xlabel('�ngulo (/Theta_0)');
ylabel('K(/Theta_0)');
legend('Fejer','Clenshaw-Curtis');
hold off;
diferencia = abs(K_theta_fej - K_theta_cc);
figure;
semilogy(valorsTheta, diferencia);
title('Diferencia entre Clenshaw-Curtis i F�jer');
xlabel('�ngulo (/Theta_0)');
ylabel('abs(K_/Theta_f_e_j - K_/Theta_c_c');

%{
Observo que per a valors de theta compresos entre 0 i 2 la diferencia
entre el resultat de la integral amb F�jer i Clenshaw Curtis �s
gaireb� zero (s'arriba a precisi� de m�quina). 

Tanmateix, a mesura que la theta va augmentant el seu valor i 
acostant-se a pi, la cuadratura de F�jer comen�a a distanciar-se 
de la de Clenshaw-Curtis de manera que, quan theta �s gaireb� pi, 
la difer�ncia entre les dues cuadratures �s m�xima.

Ara ens disposem a comprobar l'error de la cuadratura de F�jer quan ens
acostem a pi segons el nombre de punts utilitzats:
A l'eix X hi posarem les m per les quals evaluem i al Y hi posarem el
error maxim dels ultims punts m�s propers a pi ([pi-pi/2:10^-2:pi-10^-15]).
%}
thetas = [0:0.99*pi/100:pi];
propersAPi = sin(thetas/2);
Ms = [0:6:450];
anterior = ones(1, length(propersAPi));
actual = zeros(1, length(propersAPi));
error = zeros(1, length(propersAPi));
Ys = [];
mNes = [];
%{
for theta = thetas
   m = 2;
   error = 1;
   actual = 1;
   anterior = 1;
   while error > 10^-10
       k = sin(theta/2);
       K = @(phi)(1./((1-(k.^2).*(sin(phi)).^2)).^(1/2));
       %Aqui sembla que no hi arriba mai.
       actual = fejer(a, b, m, K);
       error = abs(actual - anterior);
       anterior = actual;
       m = m + 2;
   end
    mNes = [mNes m];
end

figure();
plot(thetas, mNes);
title('m necessari per obtenir error 10^-10');
xlabel('Angle (/Theta)');
ylabel('m necessari per obtenir error 10^-10');


for m = Ms
    comptador = 0;
    for k = propersAPi
        comptador = comptador + 1;
        K = @(phi)(1./((1-(k.^2).*(sin(phi)).^2)).^(1/2));
        actual(comptador) = fejer(a, b, m, K);
    end
    error = abs(actual - anterior);
    anterior = actual;
    Ys = [Ys max(error)];
end
figure();
plot(Ms, Ys);
title('Error absolut maxim prop de pi respecte m');
xlabel('m');
ylabel('Error absolut maxim prop de pi');
%}

% Tanmateix si fem salts de 10 en 10 trobem un resultat totalment diferent:
thetas = [pi-pi/2:10^-2:pi-10^-15];
propersAPi = sin(thetas/2);
Ms = [50:6:450];
anterior = ones(1, length(propersAPi));
actual = zeros(1, length(propersAPi));
error = zeros(1, length(propersAPi));
Ys = [];

for m = Ms
    comptador = 0;
    for k = propersAPi
        comptador = comptador + 1;
        K = @(phi)(1./((1-(k.^2).*(sin(phi)).^2)).^(1/2));
        actual(comptador) = fejer(a, b, m, K);
    end
    error = abs(actual - anterior);
    anterior = actual;
    Ys = [Ys max(error)];
end
figure();
plot(Ms, Ys);
title('Error absolut maxim prop de pi respecte m');
xlabel('m');
ylabel('Error absolut maxim prop de pi');


%Quan theta = pi:
theta = pi;
k = sin(theta/2);
K = @(phi)(1./((1-(k.^2).*(sin(phi)).^2)).^(1/2));
K_theta_fejPI = fejer(a, b, m, K);
K_theta_ccPI = cuadratura_cc(a, b, m, K);
K_theta_fejPI
K_theta_ccPI
%{
Tal com ens diu l'enunciat, quan theta = pi la integral �s
impr�pia. En interpretar el sentit f�sic que aix� t�, dedu�m que 
aquesta integral NO �s convergent, ja que si l'angle des del
qual �s deixat anar (amb velocitat nula) el p�ndol �s pi, �s a dir,
que el p�ndol es troba en posici� vertical, aleshores el seu periode no
pot donar cap valor concret, donat que el p�ndol es quedar� quiet
en aquesta posici� i aix� mai arribar� a oscil�lar.

Observem que la quadratura que ens aporta un resultat l�gic �s 
la de Clenshaw_Curtis, que proporciona -Inf com a resultat.
La cuadratura de F�jer, en canvi, ens d�na un valor exacte i aix�
erroni. Tal com hem pogut observar en els gr�fics anteriors, aquesta
quadratura comen�a a funcionar malament a partir de theta = 2 i
arriba al seu error m�xim quan theta~pi, que �s quan la integral es fa
impr�pia.
%}