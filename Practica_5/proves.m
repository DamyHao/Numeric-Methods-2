%% PRACTICA_15_Casas_Mercade


%% Section A)
% For the temperatures T = {0.99, 0.98, 0.97, . . . , 0.85}, use Newtons
% method to compute the coordinates vl(T) and vg(T), along with the
% corresponding pressure. When changing T, use the previous result as
% initial guess. For T = 0.99, use v(0) = 0.8 and v(0) = 1.2.
format long
close all;
clear all;

Ts = 0.99:-0.01:0.99; %Vector of temperatures

v0 = [0.8; 1.2]; %Initial guess for T=0.99

v = [v0]; %vector of volumes
V = []; %matrice that will contain the final solution for each temperature

xarxa = 0.3:0.005:5; %the volume points where the fucntion will be evaluated in order to plot the isotherm lines 

colors = jet(length(Ts)); % Rainbow colors for the plot 
colors = fliplr(colors); %We fliplr in order to have the red colors for highest temperatures and blue for the lowest
ii = 0;
figure(1);

for T = Ts
    %The functionP5 have the two equations that have to be solved
    funcioP5 = @(x)([log((3 * x(2) - 1) / (3 * x(1) - 1)) + 9 / (4 * T) * (1 / x(2) - 1 / x(1)) - 1 / (3 * x(2) - 1) + 1 / (3 * x(1) - 1), (8 * T) / 3 * (1 / (3 * x(2) - 1) - 1 / (3 * x(1) - 1)) - 1 / x(2)^2 + 1 / x(1)^2]);
    vanDerWaals = @(vol)((8 * T ./ (3 .* vol - 1)) - (3 ./ (vol.^2)));
    ii = ii + 1;
    plot(xarxa, vanDerWaals(xarxa)); %here we represent an isotherm line
    text(1.1, vanDerWaals(1.1), ['T = ', num2str(T)], 'FontSize', 6); %It shows which isotherm is showing the plot
    axis([0.5, 3.5, 0, 1.5]);
    hold on;

    [XK, resd, it] = newtonn(v, 10^-8, 100, funcioP5); %we call newtonn to solve funcioP5
    v = XK(:, end); %We take the result to be the new initial guess for the next temperature
    V = [V, v]; %we save this result in the matrice

    plot(v(1), vanDerWaals(v(1)), '-*', 'Color', colors(ii, :)); %plot of the point (vl, pl)
    hold on

    plot(v(2), vanDerWaals(v(2)), '-*', 'Color', colors(ii, :)); %plot of the point (vg,pg)
    hold on

end

set(get(gca, 'XLabel'), 'String', 'v');  %description of each axis
set(get(gca, 'YLabel'), 'String', 'P');
hold off