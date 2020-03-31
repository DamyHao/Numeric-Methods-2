clear all
close all
clc

%addpath('../Practica_5')

%% Section A)

determinants = [];
alphas = 0:0.01:3;

alphaZeros = [];
posicio=[];

figure(1)
i=1;

for alpha=alphas
    
    f=@(phi)([tan(phi(1))-alpha*(2*sin(phi(1)) +sin(phi(2))) ; tan(phi(2)) - 2*alpha*(sin(phi(1)) + sin(phi(2)))]);
    
    phi=[0,0];
    
    j = jaco(f,phi);
    
    determinants = [determinants, det(j)];
    
    if abs(det(j)) < 0.01
        alphaZeros = [alphaZeros, alpha];
        posicio=[posicio i];
    end
    i=i+1;
end

Det0=[determinants(posicio(1)), determinants(posicio(2))];
plot(alphas, determinants,'LineWidth',2)
hold on
title('Det(J(0,0)) as a function of alpha')
xlabel('Alpha')
ylabel('Det(J(0,0))')
plot(alphaZeros,Det0,'*')
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

alphas = 0:0.1:2;
% Dominis dels angles
dom1 = [0, pi/2];
dom2 = [-pi/2, pi/2];
factor = [dom1(2); (dom2(1)-dom2(2))];
a = [dom1(1); dom2(1)];
aleatoryTimes = 1:20;

sol = [];
alphasol=[];
figure(2)
for alpha = alphas
    f=@(phi)([(tan(phi(1))-alpha*(2*sin(phi(1)) +sin(phi(2)))) , (tan(phi(2)) - 2*alpha*(sin(phi(1)) + sin(phi(2))))]);
    
    
    for i = aleatoryTimes
        aleatory = rand(2,1);
        phi0=aleatory.*factor-a;
        %x*(a-b) - a
        
        [XK, resd, it] = newtonn(phi0, 1e-6, 100, f);
        % Comprobar que estigui dintre el domini
        if XK(1, end) > dom1(1) && XK(1,end) < dom1(2) && XK(2,end) > dom2(1) && XK(2, end) < dom2(2)
            % Comprobar que no sigui 0:
            if XK(1, end) ~= 0 XK(2,end) ~= 0
                sol=[sol, XK(:,end)];
                alphasol=[alphasol, alpha];
            end

         
        end
         
         
    end
    
    
end

subplot(2,1,1)
plot(alphasol, sol(1,:),'o','Color','blue')
axis([0 2 -0.8 1.6])
title('Phi1 as a function of alpha')
xlabel('Alpha')
ylabel('Phi1')
subplot(2,1,2)
plot(alphasol, sol(2,:),'o','Color','y');
axis([0 2 -0.8 1.6])
title('Phi2 as a function of alpha')
xlabel('Alpha')
ylabel('Phi2')

figure (3)
plot(alphasol, sol(2,:),'o','Color','y');
hold on
plot(alphasol, sol(1,:),'o','Color','blue')
axis([0 2 -0.8 1.6])
title('Angles of rotation as a function of alpha')
legend('Phi2','$\hat{\psi}$','Location','southwest','Interpreter','latex')
xlabel ('Alpha')
ylabel('Phi')
hold off
%<<mates.jpg>> 


%img = imread('mates.jpg');
%image(img);

%% Section c)

funAlpha =@(y)([tan(y(1))-y(3)*(2*sin(y(1)) +sin(y(2))) ; tan(y(2)) - 2*y(3)*(sin(y(1)) + sin(y(2)))]);
%Agafem les primeres solucions de l'esquerra:
y0 = [sol(:, 1); alphasol(1)];

y1 = [sol(:, 2); alphasol(2)];
y = y1; % Nomes per la entrada al while
s = 1;
tol = 1e-6;
itmax = 100;
Y=[];

while s > 0 && y(3) < 2.1 && y(3) > 0
    [y, iconv] = continuationStep(funAlpha, y0, y1, s, tol, itmax);
 
    if iconv == 1 % No hem aconseguit solució
        s = s - 0.1;
    else
        y0 = y1;
        y1 = y;
        Y=[Y , y]; %solucions
    end
    y(3)
end
size(Y)

figure(4);
plot(Y(3, :), Y(1,:));