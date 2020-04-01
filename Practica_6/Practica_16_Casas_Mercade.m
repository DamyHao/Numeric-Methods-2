clear all
close all
clc

%addpath('../Practica_5')

%% Section A)

determinants = [];
alphas = 0:0.01:3;

alphaZeros = [];
posicio = [];

figure(1)
i = 1;

for alpha = alphas

    f = @(phi)([tan(phi(1)) - alpha * (2 * sin(phi(1)) + sin(phi(2))); tan(phi(2)) - 2 * alpha * (sin(phi(1)) + sin(phi(2)))]);

    phi = [0, 0];

    j = jaco(f, phi);

    determinants = [determinants, det(j)];

    if abs(det(j)) < 0.01
        alphaZeros = [alphaZeros, alpha];
        posicio = [posicio i];
    end

    i = i + 1;
end

Det0 = [determinants(posicio(1)), determinants(posicio(2))];
plot(alphas, determinants, 'LineWidth', 2)
hold on
title('Det(J(0,0)) as a function of alpha')
xlabel('Alpha')
ylabel('Det(J(0,0))')
plot(alphaZeros, Det0, '*')
% The implicit function theorem (imft) states that as long as the jacobian
% is non-singular (det non zero) the system will define phi(1) and phi(2)
% as a unique functions of aplha, so we'll have a unique map between
% the solutions and alphas
%When the determinant is zero the uniqueness will be lost locally nearby
%the aplha points which make the determinant 0 and new branches of
%solution may emerge.
disp('The alpha values that make zero the determinant are more less:')
disp(alphaZeros)

%% Section B)
%addpath('..Practica_5')
x = 0.01;
alphas = 0:x:2;
% Dominis dels angles
dom1 = [0, pi / 2];
dom2 = [-pi / 2, pi / 2];
factor = [dom1(2); (dom2(1) - dom2(2))];
a = [dom1(1); dom2(1)];
aleatoryTimes = 1:20;

sol = [];
alphasol = [];
figure(2)

for alpha = alphas
    f = @(phi)([(tan(phi(1)) - alpha * (2 * sin(phi(1)) + sin(phi(2)))), (tan(phi(2)) - 2 * alpha * (sin(phi(1)) + sin(phi(2))))]);

    for i = aleatoryTimes
        aleatory = rand(2, 1);
        phi0 = aleatory .* factor - a;
        %Formula per canviar de escala i moure:
        %x*(a-b) - a

        [XK, resd, it] = newtonn(phi0, 1e-6, 100, f);
        % Comprobar que estigui dintre el domini
        if XK(1, end) > dom1(1) && XK(1, end) < dom1(2) && XK(2, end) > dom2(1) && XK(2, end) < dom2(2)
            sol = [sol, XK(:, end)];
            alphasol = [alphasol, alpha];
        end

    end

end

subplot(2, 1, 1)
plot(alphasol, sol(1, :), 'o', 'Color', 'blue')
axis([0 2 -0.8 1.6])
title('Phi1 as a function of alpha')
xlabel('Alpha')
ylabel('Phi1')
subplot(2, 1, 2)
plot(alphasol, sol(2, :), 'o', 'Color', 'y');
axis([0 2 -0.8 1.6])
title('Phi2 as a function of alpha')
xlabel('Alpha')
ylabel('Phi2')

figure(3)
plot(alphasol, sol(2, :), 'o', 'Color', 'y');
hold on
plot(alphasol, sol(1, :), 'o', 'Color', 'blue')
axis([0 2 -0.8 1.6])
title('Angles of rotation as a function of alpha')
legend('Phi2', 'Phi1', 'Location', 'southwest', 'Interpreter', 'latex')
xlabel('Alpha')
ylabel('Phi')
hold off


%% Section c)


%Agafem les primeres solucions de la drieta:

%As it has been seen analitically there are two alpha's that make zero the
%jocobian determinant, this alphas are 0.293 and 1.707, so the in order to
%obtand all the solutions with the secant continuation step we will take 3
%difrent y0 and y1 at the right plot of the anterious exercise
funAlpha = @(y)([tan(y(1)) - y(3) * (2 * sin(y(1)) + sin(y(2))); tan(y(2)) - 2 * y(3) * (sin(y(1)) + sin(y(2)))]);
figure(4)
epsilon = 0.01;
alphas = [2 - epsilon, 2 - 2 * epsilon];
MP = [[]];

dom1 = [0, pi / 2];
dom2 = [-pi / 2, pi / 2];

% Trobem ara dos solucions per dos alphas diferents per tirar el continuationStep:
% Ens guiem per el plot fet anteriorment per trobar les solucions diferents que volem
% Sabem que hi haura 3 solucions: La 0 (trivial), una per sobre 0.6 i una per sota. Trobarem aquestes dos ultimes
for jj = 1:2
    sol1 = 0; % Si es 0 es que falta trobarla:
    sol2 = 0;
    alpha = 2 - jj*epsilon;
    f = @(phi)([(tan(phi(1)) - alpha * (2 * sin(phi(1)) + sin(phi(2)))), (tan(phi(2)) - 2 * alpha * (sin(phi(1)) + sin(phi(2))))]);
    
    while sol1 == 0 || sol2 == 0
        aleatory = rand(2, 1);
        phi0 = aleatory .* factor - a;
        %x*(a-b) - a
    
        [XK, resd, it] = newtonn(phi0, 1e-16, 100, f);
        % Comprobar que estigui dintre el domini
        if XK(1, end) > dom1(1) && XK(1, end) < dom1(2) && XK(2, end) > dom2(1) && XK(2, end) < dom2(2)
            % I a mes que no sigui 0:
            if XK(1, end) > epsilon || XK(1, end) < (-0.001)% El blau sempre esta per sobre i nomes cal que comprobem aquest
                %Classificar si es de dalt o de sota
                if XK(1, end) > 0.6
                    MP(:, jj) = [XK(:, end); alpha];
                    sol1 = 1;
                else
                    MP(:,jj+2) = [XK(:, end); alpha];
                    sol2 = 1;
                end
            end
        end
    end
end


% Hi afegim la solució 0,0 tambe:
MP(:,5) = [0; 0; 2-epsilon];
MP(:,6) = [0; 0; 2-2*epsilon];

%Matriu dels punts on s'iniciaran els contSteps cap a la esquerra:
% Solucions:
MP

s = 1; % En un principi utilitzarem s = 1 per mantenir una separacio mes o menys similar entre les solucions.
% Aquesta separacio vindra definida per el parametre epsilon declarat anteriorment.
tol = 1e-6;
itmax = 100;
Y = [];

for it = 1:2:length(MP)% Tirarem el continuation step en les 3 branques de solucions trobades:
    y0 = MP(:, it);
    y1 = MP(:, it + 1);
    y = y1;

    while y(3) < 2 && y(3) > 0 && s > 0
        [y, iconv] = continuationStep(funAlpha, y0, y1, s, tol, itmax);

        if iconv == 1 % No hem aconseguit solució i ajustem s
            s = s - 0.1; % Si la s arriba a 0 desistirem i no buscarem mes solucions
        else
            y0 = y1;
            y1 = y;
            % Nomes les guardem si estan dintre el domini
            if y(1, end) > dom1(1) && y(1, end) < dom1(2) && y(2, end) > dom2(1) && y(2, end) < dom2(2)
                Y = [Y, y]; %solucions
            end
        end
        plot(Y(end, :), Y(1:2, :), 'o');
        hold on
    end
end

hold off


% Dibuix del pendul en alfa = 2.14:
% Com que volem una alfa determinada utilitzarem el metode de newton amb exploracio i 
% no el del continuationStep per trobar les solucions:
figure(5)
aleatoryTimes = 1:1:1000;
alpha = 2.14;
f = @(phi)([(tan(phi(1)) - alpha * (2 * sin(phi(1)) + sin(phi(2)))), (tan(phi(2)) - 2 * alpha * (sin(phi(1)) + sin(phi(2))))]);
for i = aleatoryTimes
    aleatory = rand(2, 1);
    phi0 = aleatory .* factor - a;
    %x*(a-b) - a
    [XK, resd, it] = newtonn(phi0, 1e-16, 100, f);
    % Comprobar que estigui dintre el domini
    if XK(1, end) > dom1(1) && XK(1, end) < dom1(2) && XK(2, end) > dom2(1) && XK(2, end) < dom2(2)
            plot(alpha, XK(1, end), '*');  %phi1
            plot(alpha, XK(2,end), '*') %phi2
            hold on   
    end
end
hold off

% Canviem ara la direccio en la que tirem el continuationStep i anem cap a la dreta.
% Posarem com a maxim alfa = 3;
maxAlfa = 3;
perDibuixar = [];

figure;
% Tirarem el continuation step en les 3 branques de solucions trobades anteriorment:
for it = 1:2:length(MP)
    y1 = MP(:, it); % Per ferho cal canviar el sentit, ja que ara anirem cap a la dreta
    y0 = MP(:, it + 1);

    y = y0;
    
    sensePintar = 1;

    while y(3) < maxAlfa && y(3) > 0 && s > 0
        [y, iconv] = continuationStep(funAlpha, y0, y1, s, tol, itmax);
        if iconv == 1 % No hem aconseguit solució i ajustem s
            s = s - 0.1; % Si la s arriba a 0 desistirem i no buscarem mes solucions
        else
            y0 = y1;
            y1 = y; 
            Y = [Y, y]; %solucions
            % Pintarem el pendul solament una vegada:
            if abs(alpha-y(3)) < epsilon && sensePintar == 1
                perDibuixar = [perDibuixar, y];
                sensePintar = 0;
            end
        end
        
        plot(Y(end, :), Y(1:2, :), 'o');
        hold on
    end
    
end
hold off

% Dibuixem els pendols:
for ii = perDibuixar
    figure;
    dibuixarPendul(ii(1), ii(2), 1);
end



% Com es pot observar, les solucions continuen cap a la dreta amb normalitat. 
% No es troben noves branques ni amb el continuationStep ni amb l'exploracio newton aleatoria
% Per tant determinem que no hi ha noves branques.
% Aixo es deu a que les noves branques nomes poden neixer quan el determinant del jacobia es 0.
% Mirem com evoluciona el jacobia per alfas mes grans de dos i veiem que no torna a ser 0 (fet a l'exercici a)