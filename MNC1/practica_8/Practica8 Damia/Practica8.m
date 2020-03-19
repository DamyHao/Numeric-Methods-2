%% P8: CASAS Y JIMENEZ
%% 1) Diferencia finita local:
%% A partir de la f贸rmula de derivaci贸n aproximada:
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


%% A partir de la matriz de derivaci贸n:
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
    % No cal que toquem el primer element perq ja es 0 en els dos.
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

%% Error:
format long g;
errorsMax = [];
Ns = [10:10:1000];

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
loglog(Ns, errorsMax);
[r,m,~] = regression(log10(Ns), log10(errorsMax));

%% 2) Derivaci贸n global:
clear all;
clc;
close all;
errorsMax = [];

Ns = [4:4:32];

for n = Ns
    x = linspace(-1,3,n+1);
    % N+1 punts equiespaciats entre -1 i 3.
    fx = fun(x);
    fx = fx';
    
    D = diferenciacion(x);
    fd_aprox = D*fx;
    
    % Vector de derivades exactes:
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    fd_aprox = fd_aprox';
    
%     if n == 4 || n == 32
%         figure;
%         plot(x, fd);
%         hold on;
%         plot(x,fd_aprox);
%         hold off;
%     end
    error = max(abs(fd_aprox-fd));
    errorsMax = [errorsMax error];
end

figure;
loglog(Ns, errorsMax);

% En aquesta figura 閟 dificil veure Runge, tanmateix el podem observar en
% el rebot de l'error que es produeix a partir de n = 28. Si fem el mateix
% proc閟 pero arribant a n m閟 altes (per exemple 128) podem veure el
% fenomen de Runge m閟 clarament:

clear all;
%clc;
close all;
errorsMax = [];

Ns = [4:4:128];

for n = Ns
    x = linspace(-1,3,n+1);
    % N+1 punts equiespaciats entre -1 i 3.
    fx = fun(x);
    fx = fx';
    
    D = diferenciacion(x);
    fd_aprox = D*fx;
    
    % Vector de derivades exactes:
    fd = (exp(-x)).*(sin(2.*x)).*(4.*cos(2.*x)-sin(2.*x));
    fd_aprox = fd_aprox';
    
%     if n == 4 || n == 32
%         figure;
%         plot(x, fd);
%         hold on;
%         plot(x,fd_aprox);
%         hold off;
%     end
    error = max(abs(fd_aprox-fd));
    errorsMax = [errorsMax error];
end

figure;
loglog(Ns, errorsMax);

% En aquesta gr鄁ica 閟 f郼il veure el fenomen de Runge ja que apareixen
% errors de fins a 10^25.
