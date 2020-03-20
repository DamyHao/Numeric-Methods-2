%% PRACTICA_5:NEWTON'S METHOD IN R^n
% Damia Casas & N�ria Mercad�

%% A)
% For the temperatures T = {0.99, 0.98, 0.97, . . . , 0.85}, use Newton?s
% method to compute the coordinates vl(T) and vg(T), along with the
% corresponding pressure. When changing T, use the previous result as
% initial guess. For T = 0.99, use v(0) = 0.8 and v(0) = 1.2.
format long g;
close all;
clear all;

Ts = 0.99:-0.01:0.85;

v0 = [0.8; 1.2];

v = [v0]; V = [v0];

xarxa = 0.3:0.005:5;

colors = jet(length(Ts)); % Colors arc de sant marti
ii = 0;
figure();
for T = Ts
    funcioP5 = @(x)([log((3 * x(2) - 1) / (3 * x(1) - 1)) + 9 / (4 * T) * (1 / x(2) - 1 / x(1)) - 1 / (3 * x(2) - 1) + 1 / (3 * x(1) - 1), (8 * T) / 3 * (1 / (3 * x(2) - 1) - 1 / (3 * x(1) - 1)) - 1 / x(2)^2 + 1 / x(1)^2]);
    vanDerWaals = @(vol)((8*T./(3.*vol-1))-(3./(vol.^2)));
    ii = ii + 1;
    plot(xarxa, vanDerWaals(xarxa), 'Color', colors(ii, :));
    text(1.1, vanDerWaals(1.1), ['T = ', num2str(T)], 'FontSize', 6); %Text a cada isotermica
    axis([0.5, 3.5 , 0, 1.5]);
    hold on;


    [XK, resd, it] = newtonn(v, 10^-8, 100, funcioP5);
    v = XK(:, end);
    V = [V, v];

    plot(v(1), vanDerWaals(v(1)), '-*', 'Color', colors(ii, :));
    hold on

    plot(v(2), vanDerWaals(v(2)), '-*', 'Color', colors(ii, :));
    hold on



end

set(get(gca, 'XLabel'), 'String', 'v');
set(get(gca, 'YLabel'), 'String', 'P');
hold off

% b)

