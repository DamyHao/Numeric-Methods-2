%% P11: CASAS  Y JIMENEZ
%% 1) Integración de K_theta para theta = pi/6, pi/4, pi/3, pi/3 i 3*pi/4:
clear all;
format long g;
clc;

%{
Para resolver este apartado utilizaremos tres cuadraturas distintas, 
lo cual nos permitirá dar mayor validez a nuestro resultado y
contrastar y comparar su funcionamiento para integrar:

K = @(phi)(1./((1-(k.^2).*(sin(phi)))).^(1/2));

Se resuelve la integral mediante:

- Trapezoidal Compuesta: Dado que funciona especialmente bien para
funciones periódicas como ésta (aparece un seno).

- Cuadratura de Clenshaw - Curtis: Dado que se fundamenta en los
polinomios de Chebyshev, más estables numéricamente, y que permite así
obtener un buen resultado para la mayoría de casos.

- Cuadratura de Féjer: Dado que se corresponde con la cuadratura de
Clenshaw - Curtis abierta (sin incluir los extremos), y que aplicada sobre 
integrales que presentan alguna singularidad no es demasiado exacta, 
pero nos puede proporcionar una primera aproximación al valor de la 
integral que nos sirva de referencia a la hora de aplicar otros métodos.

Las funciones que calculan el valor de la integral mediante los
métodos anteriormente presentados son, respectivamente, las siguientes:

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
%   La funcio entregada fun ha de ser capaÃ§ de tractar vectors si se li
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
        % n serÃ  sempre par per aixo podem dividir per 2.
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
%   La funcio entregada fun ha de ser capaÃ§ de tractar vectors si se li
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

%Límites de integración:
a = 0;
b = pi/2;

%Requerimiento del problema: 10 cifras significativas
tol = 10^-10;

valorsTheta = [pi/6, pi/4, pi/3, pi/3, 3*pi/4];
for theta = valorsTheta
    %La funció a integrar s'actualitza segons el valor de theta:
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
    while error > tol %FÉJER
        I_fej = fejer(a, b, m3, K);
        error = abs(actual - I_fej);
        actual = I_fej;
        m3 = m3+10;
    end
%Y proporciono los valores que voy obteniendo para cada valor de theta:
display = sprintf('Para theta = %d, el valor de K, según el método utilizado, es: I_trap_c = %d, I_cc = %d, I_fej = %d',theta,I_trap_c,I_cc,I_fej);
disp(display);
end
%% 2) Aproximación del periodo para pequeñas oscilaciones. Ángulo a partir del cual el error de la aproximación supera 5%:
g = 9.81;
l = 1;
%Aproximación del periodo del péndulo para pequeñas oscilaciones:
T_posc = 2*pi*((l/g)^(1/2));

%{
Primera manera de hacerlo:
Creo un bucle que vaya aumentando el valor de theta progresivamente en
función de una precisión fijada, e impongo que ese bucle 
finalice cuando la diferencia entre el valor del periodo obtenido
mediante la fórmula de pequeñas oscilaciones y el obtenido mediante
la fórmula real sea superior al 5% del valor real del periodo (es decir,
del valor obtenido con la fórmula que se debe integrar).
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
%% 3) Representación de K_theta en función del ángulo:
%{
En este apartado deberemos calcular y representar K_theta desde theta=0
hasta valores muy proximos al ángulo pi, para el cual prevemos que
la integral presentará una singularidad de tipo asintótico.

Utilizaremos únicamente los métodos de Féjer, que nos proporcionará
una primera aproximación a la integral, y el de Clenshaw-Curtis, que
proporcionará un resultado previsiblemente mejor. 

Observaremos también que el método de Féjer no funciona bien en
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
title('Valor de la integral según el ángulo');
xlabel('Ángulo (theta_0)');
ylabel('K(theta_0)');
hold off;
diferencia = abs(K_theta_fej-K_theta_cc);
figure;
semilogy(valorsTheta, diferencia);
title('Diferencia entre Clenshaw-Curtis i Féjer');
xlabel('Ángulo (theta_0)');
ylabel('abs(K_theta_fej - K_theta_cc)');

%{
Observo que per a valors de theta compresos entre 0 i 2 la diferencia
entre el resultat de la integral amb Féjer i Clenshaw Curtis és
gairebé zero (s'arriba a precisió de màquina). 

Tanmateix, a mesura que la theta va augmentant el seu valor i 
acostant-se a pi, la cuadratura de Féjer comença a distanciar-se 
de la de Clenshaw-Curtis de manera que, quan theta és gairebé pi, 
la diferència entre les dues cuadratures és màxima.
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
Tal com ens diu l'enunciat, quan theta = pi la integral és
impròpia. En interpretar el sentit físic que això té, deduïm que 
aquesta integral no ha de ser convergent, ja que si l'angle des del
qual és deixat anar (amb velocitat nula) el pèndol és pi, és a dir,
que el pèndol es troba en posició vertical, aleshores el seu periode no
pot donar cap valor concret, donat que el pèndol es quedarà quiet
en aquesta posició i així mai arribarà a oscil·lar.

Observem que la quadratura que ens aporta un resultat lògic és 
la de Clenshaw_Curtis, que proporciona -Inf com a resultat.
La cuadratura de Féjer, en canvi, ens dóna un valor exacte i així
erroni. Tal com hem pogut observar en els gràfics anteriors, aquesta
quadratura comença a funcionar malament a partir de theta = 2 i
arriba al seu error màxim quan theta~pi, que és quan la integral es fa
impròpia.
%}