%% P4: Casas y Jimenez
%% a) Plantear la ecuación de la intersección del montículo y la parábola, y representar gráficamente las tres funciones:
clear all;
g = 9.81; a = 0.15; b = 0.04; %son las constantes.
x = [0:1/2:200]; Vo = 37; Ao = 67.5; %son las condiciones iniciales del problema.

Vxo = Vo*cosd(Ao);
Vyo = Vo*sind(Ao);

Ypar = ((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2)); %es la ecuación de la parábola de la trayectoria.
Yper = a*(x.^2).*exp((-b).*x); %es la ecuación del perfil del montículo.
dif = Ypar - Yper; %es la diferencia de las dos, que muestra sus puntos de corte.

figure(1); %es la representación gráfica de las tres ecuaciones anteriores.
plot(x,Ypar,'r'); 
hold on;
plot(x,Yper,'b');
hold on;
plot(x,dif,'g');
axis([0 200 0 60]);
title('Intersección');
xlabel('Desplazamiento(m)');
ylabel('Altura(m)');
hold off;

%% b) Determinar el punto de impacto usando el método de Newton:
%{ 
Para determinar el punto de impacto (Ximp) de la pelota con el montículo,
tendremos que hacer de la ecuación diferencia anterior una función.m de
Matlab, para así poder estudiar su convergencia mediante el método de
Newton, el cual tendremos localizado en otro archivo.
 %}

diferencia = @(x)((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2)) - (a*(x.^2).*exp((-b).*x));
[Ximp,it] = newton(60,120,10^-10,20,diferencia);
%{  
Considero ini=60 i fin=120 dado que, según la gráfica del apartado a)
(figure(1)), el punto de corte que me interesa se sitúa entre estos dos
valores, tal y como se puede ver a simple vista.
 %} 

%{  
La función Newton utilizada en este apartado es la siguiente:

function [Ximp,it] = newton(ini,fin,tol,itmax,fun)
it=0; Ximp=[(ini+fin),(ini+fin)/2,(ini+fin)/3]; tolk=1;
  if fun(ini)*fun(fin)<0 
      while tolk>tol && it<itmax
          Ximp = [Ximp Ximp(end)-(fun(Ximp(end))/der(Ximp(end),fun))];
          it = it+1;
          tolk = abs(Ximp(end)-Ximp(end-1));
      end
  else
      disp('Les dues imatges tenen el mateix signe.');
  end
end
 %}
%{ 
La derivada que utiliza Newton en su 5ª línea se halla en un archivo der.m
de la forma:

function d = der(Xo,f)
  h= 10^-10;
  d = (f(Xo+h)-f(Xo))/h;
end
 %}

%% c) Representar gráficamente el alcance y el error en función de la iteración:
[xk, err, it] = newton_err(60,120,10^-10,20,diferencia);

eixX = [0:1:it];
figure(2);
plot(eixX, xk, '-oy');
title('Alcance');
xlabel('Iteraciones');
ylabel('Alcance (m)');

figure(3);
plot(eixX, err, '-ok');
title('Error');
xlabel('Iteraciones');
ylabel('Error (m)');
axis([2.5 7.5 (-0.5)*10^-3 2*10^-2]);

%{  
La función Newton utilizada en este apartado es la siguiente:

function [xk,err,it] = newton_err(a,b,tol,itmax,fun)
it=0; iterr=1; tolk=1; xk = [(a+b)/2];

    if fun(a)*fun(b)<0
        while (it<itmax && tolk>tol)
            xk = [xk xk(end)-(fun(xk(end))/der(xk(end),fun))];
            tolk = abs(xk(end-1)-xk(end));
            it = it+1;
        end
        err = [];
        while iterr < it + 2
            err = [err abs(xk(iterr) - xk(end))];
            iterr = iterr + 1;
        end
                
    else
        disp('Les dues imatges tenen el mateix signe.');
    end
end
%} 

%% d) Orden de convergencia del Newton:

[Ximp,it] = newton(60,120,10^-10,20,diferencia);
%{
Vuelvo a calcular Ximp entre los mismos puntos (no hace falta, pero lo hago para recordar que sería necesario 
si se ejecutase este apartado de manera aislada). 
%}

figure(4)
EK = abs(Ximp(1:end-1)-Ximp(end)); %defino las épsilon sub k.
Xek = log10(EK(1:end-1)); %establezco Xek como las épsilon en k. 
Yek = log10(EK(2:end)); %y Yek como las épsilon en k+1.
plot(Xek,Yek,'-ob');
title('Orden de convergencia');
xlabel('log E(k)');
ylabel ('log E(k+1)');
hold on
n=(-10:0.1:2); %y calculo ahora la pendiente de la función.
m=(1/1.8).*n;
plot(m,n,'r');

%{ 
La función Newton utilizada en este apartado es la misma que en el apartado
b), es decir, la que se limita a calcular el punto de convergencia.
%}

%% e) Variamos Ao en el rango [5º,80º]:
%% i) Representar gráficamente Ximp y it en función de Ao:
%{ 
La función Newton utilizada en el apartado e) es newtonD y es la siguiente:
function [xk,it,p] = newtonD(a,b, tol, itmax, fun)
it=0; 
tolk=1;
xk = [(a+b)/2];

    if fun(a)*fun(b)<0
        while (it<itmax && tolk>tol)
            xk = [xk xk(end)-(fun(xk(end))/derD(fun,xk(end)))];
            tolk = abs(xk(end-1)-xk(end));
            it = it+1;
        end
    else
        disp('les imatges de a i b tenen signes iguals');
    end
    p = 0;
end
%}

%{ 
La función derivada también cambia i pasa a llamarse derD:
function d = derD(f,Xo);
h= 10.^-10;
f1 = f(Xo+h);
f2 = f(Xo);
d = (f1-f2)/h;
end
%}

clear all
alfaX = [5:1:80];
    g = 9.81;
    a = 0.15;
    b = 0.04;
    Vo = 37;
    impactes = [];
    iteracions = [];

for alfa = 5:1:80
   
    Vxo = Vo*cosd(alfa);
    Vyo = Vo*sind(alfa);

    func = @(x)((((Vyo/Vxo)*x)-((1/2)*g*(x^2)/(Vxo^2))) - (a*(x^2)*exp((-b)*x)));
    
    xa = 0.0001;
    xb = 10;
    % El 211 es important per detectar quan no es toquen!
    while (xb < 211) && (func(xa)*func(xb) > 0)
        xa = xb;
        xb = xb + 10;
    end
    if xb >= 210
         % disp('no es toquen');
    else
    [xk,it] = newtonD(xa, xb, 10^-5, 30, func);
    impactes(end + 1) = xk(end);
    iteracions(end + 1) = it;
    end
    
end

figure(5);
plot(alfaX, impactes, 'b');
hold off;

figure(6);
plot(alfaX, iteracions, 'g');
hold off;

%% ii) Angulo inicial mínimo:
% El cim serà el primer punt en que la derivada sigui 0 o negativa.
% Es més eficient fer servir els resultats de l'exercici anterior:
derivadaMonti = 1;
anglesPossibles = 80-5+1;
compt = 0;
alfaX = [5:1:80];

a = 0.15;
b = 0.04;
fMonti = @(x)(a*(x^2)*exp((-b)*x));

while derivadaMonti > 0 && compt <= anglesPossibles
    compt = compt + 1;
    derivadaMonti = derD(fMonti, impactes(compt));
end

angleMinim1 = alfaX(compt);

%% iii) Irregularidad:
%{
Prop del angle mínims s'observa una irregularitat en la grafica ja que fa
un salt en el punt d'impacte important. Passa de caure en el punt 28.92 amb
un angle de 61º a el punt 97.46 si la llencem amb un angle de 62º. Per tant
només amb l'ajuda de la gràfica ja podriem haber sapigut quan la pilota
travesava la muntanya.
%}

%{
No ho podem demostrar ja que no hem realitzat un mètode iteratiu com a tal
per trobar-lo
%}

%{
El mal condicionament del problema provoca que les gráfiques siguin molt
més "abruptes". Es a dir que tenen una pendent molt alta. Aixo provoca que
el metode de Newton convergeixi molt més lent. Cal recordar que el Newton
funcionava de manera que traçava una tangent a la gràfica i buscava on
queia la tangent per afegir una nova Xk a la llista. Si la tangent té una
pendent de valor absolut molt alt (igual que la gràfica), avançarà molt poc
en la direcció de l'eix. Això provoca que convergeixi més lent que si la
grafica esta ben condicionada; es a dir que els increments de x i de y són
similars i no existeixen per tant pendents elevades.
%}
