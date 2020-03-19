clear all;
%% Exercici a);
g = 9.81;
a = 0.15;
b = 0.04;

x = [0:1/100:200];
Vo = 37;
A = 67.5;

Vxo = Vo*cosd(A);
Vyo = Vo*sind(A);

Ypar = ((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2));
Yper = a*(x.^2).*exp((-b).*x);
D = Ypar - Yper;

figure(1);
plot(x,Ypar,'r');
hold on;
plot(x,Yper,'b');
hold on;
plot(x,D,'g');
axis([0 200 0 60]);
hold off;

%% Exercici b)
func = @(x)((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2)) - (a*(x.^2).*exp((-b).*x));
[xk,it,p] = newton(70, 80, 10^-3, 200, func);

%% Exercici c)
[xk, err, it, p] = newtonRepr(70, 80, 10^-5, 200, func);

eixX = [0:1:it];
figure(2);
plot(eixX, xk, 'y');
title('Alcance');
xlabel('Iteracions');
ylabel('Alcance (m)');

figure(3);
plot(eixX, err, 'r');
title('Error');
xlabel('Iteracions');
ylabel('Error (m)');
%axis([0 5 0 10^-8]);

%% Exercici d)
clear all
func = @(x)((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2)) - (a*(x.^2).*exp((-b).*x));



%% Exercici e)
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

%% ii) 
% El cim serà el primer punt en que la derivada sigui 0 o negativa.

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

angleMinim = alfaX(compt);
% disp(angleMinim);

%% iii) 
%{
Prop del angle mínims s'observa una irregularitat en la grafica ja que fa
un salt en el punt d'impacte important. Passa de caure en el punt 28.92 amb
un angle de 61º a el punt 97.46 si la llencem amb un angle de 62º. Per tant
només amb l'ajuda de la gràfica ja podriem haber sapigut quan la pilota
travesava la muntanya.
%}
