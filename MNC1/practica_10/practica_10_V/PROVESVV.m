%% P10: CASAS Y JIMENEZ
%% 1) Potencial en el eje:
%{
Se trata de aproximar el valor del potencial gravitatorio en el 
eje de un anillo a partir de la resolución de una integral mediante
la cuadratura de Clenshaw-Curtis.
%}

%{
La función que aproxima la integral mediante Clenshaw-Curtis es 
la siguiente:

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
%}

%{
Para que, tal como nos dice el enunciado, la función anterior nos
aproxime la integral para cualquier valor (y,z)~=(1,0) debemos llamar 
a la función sustituyendo las variables por sus valores

if y~=1 && z~=0
dV = @(th,dens,radi,x,y,z)
else
disp('L'enunciat ens imposa que (y,z)~=(1,0)');
end

y después, en una malla de nodos de Chebyshev para N=20 que nos hayamos
construido según a y b (variables de integración)

j = [0:1:N];
xcheb = cos(j.*pi./N);
a=0; b=2*pi;
xk = a + ((b-a)./2).*(xcheb+1);

evaluar la función en los nodos (fx) 

fx = dV(xk,1,1,0,0,0);

y llamar a la función de Clenshaw-Curtis comentada anteriormente, 
que nos devolverá la integral aproximada:

V = clenshawcurtis_p10(a, b, N, fx)
%}

%Particularizamos el caso comentado anteriormente para y=0:
clear all;
close all;
format long g;

%Calculo los nodos de Chebyshev para aplicar Clenshaw-Curtis:
N=20;
j = [0:1:N];
xcheb = cos(j.*pi./N);
a=0; b=2*pi;
xk = a + ((b-a)./2).*(xcheb+1);

%La función que integraré es la siguiente:
dV = @(th,dens,radi,x,y,z)(-(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z^2)));

%Y éste es su valor en los nodos, con las demás variables definidas:
th = xk;
dens = 1; radi = 1;
x = 0; y = 0; z = 0;
fx = dV(th,dens,radi,x,y,z);

%La integral aproximada la obtendremos mediante:
V = clenshawcurtis_p10(a, b, N, fx)

%{
Sabemos que su resultado exacto (que se puede obtener analíticamente) es
V = -2*pi, y así estudiamos el error (o exactitud) de nuestra
aproximación restándola del valor exacto. Dado que el error obtenido es
de orden de 10^-15, podemos comprobar que nuestra aproximación es
muy precisa.
%}
format short g;
error = abs((-2*pi)-V)
%% 2) Potencial en el plano:
format long g;
Ys = [0.25 0.75 1.5];
figure;
for y = Ys
    V=[];
    Zs = [-5:0.1:5];
    for z = Zs
        %{
        Tenemos que N=20, dens=1, radi=1, a=0, b=2*pi y x=0 en todo
        el problema, por lo que evito volverlos a escribir. Además, como
        N no varía, los nodos (xk) tampoco, y así no los debo recalcular.
        %}
        th = xk;
        fx = dV(th,dens,radi,x,y,z);
        V = [V clenshawcurtis_p10(a, b, N, fx)];
    end
    plot(Zs,V);
    hold on;
end
title('Potencial en el plano y - z');
xlabel('Eje vertical (z)');
ylabel('Valor aproximado del potencial (V)');
legend('y = 0.25','y = 0.75','y = 1.5');
hold off;

%% 3) Curvas equipotenciales sobre el plano:
xo = 0;
%Los puntos de partida:
yo = 0.25; 
zo=0;

h=0.01;

Ns = 1:3000;
y=zeros(1,length(Ns)); y(1)=yo;
z=zeros(1,length(Ns)); z(1)=zo;
for n=Ns
[y(n+1),z(n+1)] = nextPoint(y(n),z(n),xk);
end
figure;
title('Prueba de la funcion nextpoint');
plot(y,z);
%{
La función que calcula el siguiente punto equipotencial es:
%}

% Repitiendo el procedimiento para otros puntos de partida:
figure;
for yo=0:0.25:1.25
zo=0;
h=0.01;

Ns = 1:3000;
y=zeros(1,length(Ns)); y(1)=yo;
z=zeros(1,length(Ns)); z(1)=zo;
for n=Ns
[y(n+1),z(n+1)] = nextPoint(y(n),z(n),xk);
end
plot(y,z);
hold on;
end
hold off;