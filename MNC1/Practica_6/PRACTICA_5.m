%% P5: CASAS Y JIMENEZ
%% a) Determinar matriu de polinomis cardinals de Lagrange:
clear all;
n = 10; 
ncont= [0:1:n];
x = -1+(2*ncont)/n; %considero 11 nodes equiespaiats.
z = [-1:1/300:1];
P = MLagrange(z, x); 

%{
La funci� que permet generar la matriu dels polinomis �s la seg�ent:

function P = MLagrange(z, x);
n = length(x); %ser� el n�mero de files.
m = length (z); %ser� el n�mero de columnes.
den = [];
num = [];
P = [;];

for i = 1:1:n  %ens genera el vector fila de la matriu
       den = x(i)-x;
       den(i) = 1; %aix� aconseguim "saltar-nos" el valor a la posici� "i".

       for k=1:1:m %ens genera el vector columna de la matriu
           num = z(k)-x;
           num(i)=1; %aix� aconseguim "saltar-nos" el valor a la posici� "i".
           P(k,i) = prod(num)/prod(den);
       end
end
%}



legend('error te�ric','error experimental','previsi� del comportament de l�error te�ric fins a n=0');
%% b) Representaci� dels polinomis de Lagrange entre nodes per a n=3,6,9:
clear all;
for n=[3,6,9]
    ncont= [0:1:n];
    x = -1+(2*ncont)/n;
    z = [-1:1/300:1];
    P = MLagrange(z, x);
    
    figure;
    for columna = 1:1:(n+1)
        plot(z, P(:,columna));
        hold on;
    end
    titol_plots_b = sprintf('Lagrange per a n = %d', n);
    title(titol_plots_b);
    ylabel('Polinomis cardinals de Lagrange');
    xlabel('Punts avaluats (nodes i z_k)');
    hold off;
end

%% c) Funci� de Lebesgue per a n=8,16,24,32:
clear all;
for n=[8,16,24,32]
    z = [-1:1/300:1];
    l = lebesgue(n);
    figure;
    semilogy(z,l); %Ajusto l'escala per poder observar el comportament de la funci�.
    titol_plots_c = sprintf('Lebesgue per a n = %d', n);
    title(titol_plots_c);
    ylabel('Valor de la funci�');
    xlabel('Punts avaluats (nodes i z_k)');
    hold on;
end
hold off;

%% d) Representaci� de l'error en la interpolaci� de f(x)=e^x per a n=4:2:60:
clear all;
format long g;
% Ja tenim feta la matriu dels polinomis cardinals de Lagrange, ara falta fer el vector amb els valors de l'exponencial.
figure;
maxError = [];
% A fora per optimitzar
z = [-1:1/300:1];
expZ = exp(z)';

for n = 4:2:60
    ncont= [0:1:n];
    x = -1+(2*ncont)/n;
    P = MLagrange(z, x);
    expX = exp(x);
    expX = expX';
    exp_transposada = P * expX;
    %{
    De manera addicional afegim la representaci� de la funci� exponencial
    interpol�lada, per verificar aix� que el nostre procediment �s el
    correcte.
    %}
    figure(8);
    plot(z,exp_transposada);
    title('Funci� exponencial interpol�lada');
    xlabel('z_k');
    ylabel('e^x');  
    hold on; 
    maxError = [maxError max(abs(expZ - exp_transposada))];
end
hold off;

figure(9);
n = [4:2:60];
semilogy(n, maxError,'r');
hold on;
llei_en = 2.^n./(n.*log(n));
semilogy(n,llei_en*10^-16, 'b');
axis([0 61 10^-16 10^5]);
hold on;
%{
Ara, per comprovar que l'extrapolaci� de E_k quan n?0 talla l'eix "y" prop de
10^-16, intentarem predir el comportament de la corba per als valors 0 fins a 4 on no est� definida.
Per fer-ho, construirem una recta entre aquests punts que tingui com a pendent
el pendent de la corba d'error te�ric en els seus dos �ltims nodes (de 4 a
6) i visualitzarem on aniria a parar si continu�s fins a n=0.
%}
punts_a_predir=0:0.1:6;
pendent_error= 0.155;
previsio_error = (10.^(pendent_error.*punts_a_predir));
%Escric la recta com a pot�ncia de 10 per contrarrestar l'efecte de l'escala logar�tmica en l'eix y.
semilogy(punts_a_predir,previsio_error*10^-16,'k');
title('Avaluaci� de l.error')
legend('error m�xim (experimental)','error te�ric','previsi� del comportament de l.error te�ric fins a n=0');
hold off;
%{ 
I, per �ltim, ampliem la gr�fica anterior per comprovar que E_k tendeix a
10^-16 quan n?0. La recta que hem construit per fer la previsi� havia de
situar-se lleugerament per damunt de la recta de l'error te�ric perqu�, si
estigues al mateix nivell, els seus valors serien lleugerament inferiors a
10^-16, la precisi� de la m�quina, i aix� no adquiriria el comportament de
recta que nosaltres desitgem.
%}
figure;
n = [4:2:60];
llei_en = 2.^n./(n.*log(n));
semilogy(n,llei_en*10^-16, 'b');
hold on;
punts_a_predir=0:0.1:6;
pendent_error= 0.155;
previsio_error = (10.^(pendent_error.*punts_a_predir));
semilogy(punts_a_predir,previsio_error*10^-16,'k');
title('Previsi� error te�ric (E_k) (vista ampliada)');
axis([0 5.5 10^-16 9*10^-16]); %canvio els eixos a una escala molt m�s reduida.
