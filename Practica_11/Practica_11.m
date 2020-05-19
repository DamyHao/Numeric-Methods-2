close all;
clc;
format long g;
%% Section A
%The analytical results are the points (0,0), (1,0), (-1,0) and the curve y=(x^3-x)/(0.9-x^2)

%% Section B
h = 1e-2;
desiredPoints = 100000;

% Using the analytical results of section A, we start the time-stepper on the following
% initial points:

figure

for guessX = -1.5:0.1:1.5
    guess = [guessX; 0];
    result = RK4(guess, h, @functionODE, desiredPoints);
    plot(result(1, :), result(2, :))
    hold on
end

hold off;

% From the results above it is clear that the system has three orbits,
% one obtained for the x points 1.5 and -1.5, another one for x=-0.5 and the last one for x=0.5.
% Meanwhile for the x points -1,0,1 there are no results obatined, because as calculated previoulsy
% in A those are equilirbrium points.

%% Section C
% Now we are looking for the stable and unstable equilibrium orbits.
% Using starting points from de previus section we will find the stable and
% unstable orbits. In ordre to find the unstable we will integrate backwards in time.

clear all;
clc
desiredPoints = 1000000;
inestableGuess1 = [-0.45; 0];
inestableGuess2 = [0.45; 0];
stableGuess = [0; 1.122];

h = 1e-2;
figure;
result = RK4(inestableGuess1, h, @functionODEBack, desiredPoints);
plot(result(1, :), result(2, :), 'LineWidth', 2, 'Color', [0.9 0.1 0.10]);
hold on;

result = RK4(inestableGuess2, h, @functionODEBack, desiredPoints);
plot(result(1, :), result(2, :), 'LineWidth', 2, 'Color', [0.9 0.1 0.10]);
hold on;

% In this case we integrate forwards to find the stable orbit:
result = RK4(stableGuess, h, @functionODE, desiredPoints);
plot(result(1, :), result(2, :), 'LineWidth', 2, 'Color', [0.4660 0.8 0.1880]);
hold on;

%For the specific intial guesses used, obatined using section B, we have found three closed trajectories, i.e periodic orbits.
% Using the forward integration we have found the stable one, and using the backwards the unstables.

%% Section D:
% We compute the Jacobian of f at the point 0,0.
% The Jacobian trace tells us if the point expands, shrinks or preserves area.

DF = jaco(@functionODE, [0; 0]);
disp('The eigenvalues for the (0,0) point are')
disp(eig(DF))
% As we expected, there is an eigenvalue that is a real positive number (unstable) and another one that is negative (stable)
% is an unstable equilibrium point

DF = jaco(@functionODE, [0.001; 0.001]);
disp('The eigenvalues for the (1,0) point are')
disp(eig(DF))

DF = jaco(@functionODE, [-0.001; -0.001]);
disp('The eigenvalues for the (-0.001; -0.001) point are')
disp(eig(DF))


guess = [0.001; 0.001];
result = RK4(guess, h, @functionODE, desiredPoints);
plot(result(1, :), result(2, :));
hold on;

result = RK4(guess, h, @functionODEBack, desiredPoints);
plot(result(1, :), result(2, :));
hold on

% We try on another point near to the origin:
guess = [-0.001; -0.001];
result = RK4(guess, h, @functionODE, desiredPoints);
plot(result(1, :), result(2, :));
hold on;

result = RK4(guess, h, @functionODEBack, desiredPoints);
plot(result(1, :), result(2, :));

legend('Forward in time, intial guess (0.001,0.001), stable orbit','Backward in time, intial guess (0.001,0.001), unstable orbit','Forward in time, intial guess (-0.001,-0.001), stable orbit','Backward in time, intial guess (-0.001,-0.001), unstable orbit', 'Location', 'best');

