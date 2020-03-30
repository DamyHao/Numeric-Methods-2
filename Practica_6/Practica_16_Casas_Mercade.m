clear all
close all
clc

addpath('../Practica_5')

%% Section A)

determinants = [];
alphas = 0:0.01:3;

alphaZeros = [];

figure(1)

for alpha=alphas
    
    f=@(phi)([tan(phi(1))-alpha*(2*sin(phi(1)) +sin(phi(2))) ; tan(phi(2)) - 2*alpha*(sin(phi(1)) + sin(phi(2)))]);
    
    phi=[0,0];
    
    j = jaco(f,phi);
    
    determinants = [determinants, det(j)];
    
    if abs(det(j)) < 0.01
        alphaZeros = [alphaZeros, alpha];
    end
    
end

plot(alphas, determinants)

% The implicit function theorem (imft) states that as long as the jacobian
% is non-singular (det non zero) the system will define phi(1) and phi(2)
% as a unique functions of aplha, so we'll have a unique map between
% the solutions and alphas
%When the determinant is zero the uniqueness will be lost locally nearby
%the aplha points which make the determinant 0 and new branches of
%solution may emerge.
disp('The alpha values that make zero the determinant are more less:')
disp(alphaZeros)

%


%% Section B)
%addpath('..Practica_5')
randomSeed = [pi/4, 0];
alphas = 0:0.1:2;
% Dominis dels angles
dom1 = [0, pi/2];
dom2 = [-pi/2, pi/2];
factor = [dom1(2); (dom2(1)-dom2(2))];
a = [dom1(1); dom2(1)];
aleatoryTimes = 1:10;

solucions = [];
figure(2)
for alpha = alphas
    f=@(phi)([(tan(phi(1))-alpha*(2*sin(phi(1)) +sin(phi(2)))) , (tan(phi(2)) - 2*alpha*(sin(phi(1)) + sin(phi(2))))]);
    
    
    for i = aleatoryTimes
        aleatory = rand(2,1);
        
        %x*(a-b) - a
        phi0 = aleatory.*factor - a;
        disp(phi0)
        [XK, resd, it] = newtonn(phi0, 1e-2, 100, f);
        
        if XK(1, end) > dom1(1) && XK(1,end) < dom1(2) && XK(2,end) > dom2(1) && XK(2, end) < dom2(2)
            plot(alpha,XK(:,end),'-*')
            hold on
            %solucions = [solucions, XK(:,end)];
        end
        
    end
    
    
end

%plot(linspace(0,2,length(solucions)), solucions);

