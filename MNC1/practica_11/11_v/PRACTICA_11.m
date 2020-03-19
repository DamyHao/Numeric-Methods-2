%% P11: CASAS  Y JIMENEZ
%% 1) Integraci�n de K_theta para theta = pi/6, pi/4, pi/3, pi/3 i 3*pi/4:
clear all;
format long g;
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
integrales que presentan alguna singularidad no es demasiado exacta, 
pero nos puede proporcionar una primera aproximaci�n al valor de la 
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
    K = @(phi)(1./((1-(k.^2).*(sin(phi)))).^(1/2));
    
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
display = sprintf('Para theta = %d, el valor de K, seg�n el m�todo utilizado, es: I_trap_c = %d, I_cc = %d, I_fej = %d',theta,I_trap_c,I_cc,I_fej);
disp(display);
end
%% 2) Aproximaci�n del periodo para peque�as oscilaciones. �ngulo a partir del cual el error de la aproximaci�n supera 5%:
g = 9.81;
l = 1;
%Aproximaci�n del periodo del p�ndulo para peque�as oscilaciones:
T_posc = 2*pi*((l/g)^(1/2));

%{
Primera manera de hacerlo:
Creo un bucle que vaya aumentando el valor de theta progresivamente en
funci�n de una precisi�n fijada, e impongo que ese bucle 
finalice cuando la diferencia entre el valor del periodo obtenido
mediante la f�rmula de peque�as oscilaciones y el obtenido mediante
la f�rmula real sea superior al 5% del valor real del periodo (es decir,
del valor obtenido con la f�rmula que se debe integrar).
%}
m = 20;
cincpercent = 1;
precisio = 0.05;
while error < cincpercent
    theta = theta + precisio;
    K_theta = fejer(a, b, m, K);
    T_real = 4*((l/g)^(1/2))*K_theta;
    error = abs(T_real - T_posc);
    cincpercent = 0.05*T_real;
end
Theta_max = theta

%segona manera de fer-ho: NO EM SURT
k = sin(theta/2);
K = @(phi)(1./((1-(k.^2).*(sin(phi)))).^(1/2));
K_theta = fejer(a, b, m, K);
T_real = 4*((l/g)^(1/2))*K_theta;
fun_theta = @(x)(abs(T_real - T_posc) - 0.05*T_real);
X = newton (2, fun_theta);
X
disp('La X no em dona el que hauria de donar');
%% 3) Representaci�n de K_theta en funci�n del �ngulo:
%{
En este apartado deberemos calcular y representar K_theta desde theta=0
hasta valores muy proximos al �ngulo pi, para el cual prevemos que
la integral presentar� una singularidad de tipo asint�tico.

Utilizaremos �nicamente los m�todos de F�jer, que nos proporcionar�
una primera aproximaci�n a la integral, y el de Clenshaw-Curtis, que
proporcionar� un resultado previsiblemente mejor. 

Observaremos tambi�n que el m�todo de F�jer no funciona bien en
acercarnos a la singularidad.
%}
valorsTheta = linspace (0, 99*pi/100, 50); %sin llegar a valer pi.
 
K_theta_cc = zeros(1,length(valorsTheta));
K_theta_fej = zeros(1,length(valorsTheta));
n = 1;
for theta = valorsTheta
    k = sin(theta/2);
    K = @(phi)(1./((1-(k.^2).*(sin(phi)))).^(1/2));
    K_theta_cc(n) = cuadratura_cc(a, b, m, K);
    K_theta_fej(n) = fejer(a, b, m, K);
    n = n + 1;
end
figure;
plot(valorsTheta, K_theta_fej);
hold on;
plot(valorsTheta, K_theta_cc);
title('Valor de la integral seg�n el �ngulo');
xlabel('�ngulo (theta_0)');
ylabel('K(theta_0)');
hold off;
diferencia = abs(K_theta_fej-K_theta_cc);
figure;
semilogy(valorsTheta, diferencia);
title('Diferencia entre Clenshaw-Curtis i F�jer');
xlabel('�ngulo (theta_0)');
ylabel('abs(K_theta_fej - K_theta_cc)');

%{
Observo que per a valors de theta compresos entre 0 i 2 la diferencia
entre el resultat de la integral amb F�jer i Clenshaw Curtis �s
gaireb� zero (s'arriba a precisi� de m�quina). 

Tanmateix, a mesura que la theta va augmentant el seu valor i 
acostant-se a pi, la cuadratura de F�jer comen�a a distanciar-se 
de la de Clenshaw-Curtis de manera que, quan theta �s gaireb� pi, 
la difer�ncia entre les dues cuadratures �s m�xima.
%}

%Quan theta = pi:
theta = pi;
k = sin(theta/2);
K = @(phi)(1./((1-(k.^2).*(sin(phi)))).^(1/2));
K_theta_fej = fejer(a, b, m, K);
K_theta_cc = cuadratura_cc(a, b, m, K);
K_theta_fej
K_theta_cc
%{
Tal com ens diu l'enunciat, quan theta = pi la integral �s
impr�pia. En interpretar el sentit f�sic que aix� t�, dedu�m que 
aquesta integral no ha de ser convergent, ja que si l'angle des del
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