%% P4: Casas y Jimenez
%% a) Plantear la ecuaci�n de la intersecci�n del mont�culo y la par�bola, y representar gr�ficamente las tres funciones:
clear all;
g = 9.81; a = 0.15; b = 0.04; %son las constantes.
x = [0:1/2:200]; Vo = 37; Ao = 67.5; %son las condiciones iniciales del problema.

Vxo = Vo*cosd(Ao);
Vyo = Vo*sind(Ao);

Ypar = ((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2)); %es la ecuaci�n de la par�bola de la trayectoria.
Yper = a*(x.^2).*exp((-b).*x); %es la ecuaci�n del perfil del mont�culo.
dif = Ypar - Yper; %es la diferencia de las dos, que muestra sus puntos de corte.

figure(1); %es la representaci�n gr�fica de las tres ecuaciones anteriores.
plot(x,Ypar,'r'); 
hold on;
plot(x,Yper,'b');
hold on;
plot(x,dif,'g');
axis([0 200 0 60]);
title('Intersecci�n');
xlabel('Desplazamiento(m)');
ylabel('Altura(m)');
hold off;

%% b) Determinar el punto de impacto usando el m�todo de Newton:
%{ 
Para determinar el punto de impacto (Ximp) de la pelota con el mont�culo,
tendremos que hacer de la ecuaci�n diferencia anterior una funci�n.m de
Matlab, para as� poder estudiar su convergencia mediante el m�todo de
Newton, el cual tendremos localizado en otro archivo.
 %}

diferencia = @(x)((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2)) - (a*(x.^2).*exp((-b).*x));
[Ximp,it] = newton(60,120,10^-10,20,diferencia);
%{  
Considero ini=60 i fin=120 dado que, seg�n la gr�fica del apartado a)
(figure(1)), el punto de corte que me interesa se sit�a entre estos dos
valores, tal y como se puede ver a simple vista.
 %} 

%{  
La funci�n Newton utilizada en este apartado es la siguiente:

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
La derivada que utiliza Newton en su 5� l�nea se halla en un archivo der.m
de la forma:

function d = der(Xo,f)
  h= 10^-10;
  d = (f(Xo+h)-f(Xo))/h;
end
 %}

%% c) Representar gr�ficamente el alcance y el error en funci�n de la iteraci�n:
[xk, err, it] = newton_err(60,120,10^-10,20,diferencia);

eixX = [0:1:it];
figure(2);
plot(eixX, xk, 'y');
title('Alcance');
xlabel('Iteraciones');
ylabel('Alcance (m)');

figure(3);
plot(eixX, err, 'r');
title('Error');
xlabel('Iteraciones');
ylabel('Error (m)');
axis([2.5 7.5 (-0.5)*10^-3 2*10^-2]);

%{  
La funci�n Newton utilizada en este apartado es la siguiente:

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
%{ 
Para determinar el orden de convergencia del m�todo de Newton en este
problema, ampliar� la funci�n Newton que he utilizado en el apartado b), la
cual me proporciona �nicamente el punto de impacto (Ximp) y el n�mero de
iteraciones, para que me proporcione tambi�n p a partir de c�lculos de
error entre iteraciones.
 %} 

[Ximp,it,Xek,Yek,p] = newton_p(60,120,10^-10,20,diferencia); %Xek es el error en n y Yek el error en n+1.
figure(4);
plot(log10(Yek),log10(Xek),'r');
title('Orden de convergencia');
xlabel('log10(Yek)');
ylabel('log10(Xek)');
axis([-8.5 0.5 0.15 1.6]);

%{ 
La funci�n Newton utilizada en este apartado es la siguiente:

function [Ximp,it,p,Xek,Yek] = newton_p(ini,fin,tol,itmax,fun) 
it=0; Ximp=[(ini+fin),(ini+fin)/2,(ini+fin)/3]; tolk=1; Xek=1; Yek=1; p=2;
  if fun(ini)*fun(fin)<0
      while tolk>tol && it<itmax
          Ximp = [Ximp Ximp(end)-(fun(Ximp(end))/der(Ximp(end),fun))];
          it = it+1;
          tolk = abs(Ximp(end)-Ximp(end-1));
          Xek = [Xek abs(fun(Ximp(end-2)))-(fun(Ximp(end-1)))];
          Yek = [Yek abs(fun(Ximp(end-1)))-(fun(Ximp(end)))];
          p = [p log10(Yek(end))/log10(Xek(end))];
      end
  else
      disp('Les dues imatges tenen el mateix signe.');
  end
end
 %} 

%% e) Variamos Ao en el rango [5�,80�]:
%% i) Representar gr�ficamente Ximp y it en funci�n de Ao:
clear all
% Cal canviar el 10 pel 5 i saber que passa en aquest interval!!
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
    [xk,it,p] = newton(xa, xb, 10^-5, 30, func);
    impactes(end + 1) = xk(end);
    iteracions(end + 1) = it;
    end
    
end

figure(4);
plot(alfaX, impactes, 'b');
hold off;

figure(5);
plot(alfaX, iteracions, 'g');
hold off;

%% ii) Angulo inicial m�nimo:
% El cim ser� el primer punt en que la derivada sigui 0 o negativa.
% Es m�s eficient fer servir els resultats de l'exercici anterior:
derivadaMonti = 1;
anglesPossibles = 80-5+1;
compt = 0;
alfaX = [5:1:80];

a = 0.15;
b = 0.04;
fMonti = @(x)(a*(x^2)*exp((-b)*x));

while derivadaMonti > 0 && compt <= anglesPossibles
    compt = compt + 1;
    derivadaMonti = der(fMonti, impactes(compt));
end

angleMinim1 = alfaX(compt);

% Tanmateix com que en el seg�ent apartat es demana analitzar la funcio per
% trobar l'angle m�nim, el trobarem mitjan�ant una funcio desde 0
% (completa). (Es podria anar m�s r�pid fent la derivada de forma
% analitica).
clear all

alfaX = [5:1:80];
    g = 9.81;
    a = 0.15;
    b = 0.04;
    Vo = 37;
    impactes = [];

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
    [xk,it,p] = newton(xa, xb, 10^-5, 30, func);
    impactes(end + 1) = xk(end);
    iteracions(end + 1) = it;
    end
    
end

%% iii) Irregularidad:
%{
Prop del angle m�nims s'observa una irregularitat en la grafica ja que fa
un salt en el punt d'impacte important. Passa de caure en el punt 28.92 amb
un angle de 61� a el punt 97.46 si la llencem amb un angle de 62�. Per tant
nom�s amb l'ajuda de la gr�fica ja podriem haber sapigut quan la pilota
travesava la muntanya.
%}

%{
No ho podem demostrar ja que no hem realitzat un m�tode iteratiu com a tal
per trobar-lo
%}

%{
El mal condicionament del problema provoca que les gr�fiques siguin molt
m�s "abruptes". Es a dir que tenen una pendent molt alta. Aixo provoca que
el metode de Newton convergeixi molt m�s lent. Cal recordar que el Newton
funcionava de manera que tra�ava una tangent a la gr�fica i buscava on
queia la tangent per afegir una nova Xk a la llista. Si la tangent t� una
pendent de valor absolut molt alt (igual que la gr�fica), avan�ar� molt poc
en la direcci� de l'eix. Aix� provoca que convergeixi m�s lent que si la
grafica esta ben condicionada; es a dir que els increments de x i de y s�n
similars i no existeixen per tant pendents elevades.
%}

