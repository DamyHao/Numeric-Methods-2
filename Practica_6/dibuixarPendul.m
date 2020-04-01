function succes = dibuixarPendul(phi1, phi2, l)
% Funcio que fa un plot de un pendul doble amb els angles indicats i la longitud de les barres indicada (les dos igual)
% a entrar els angels respecte la vertical en radiants

h = plot(0, 0, 'MarkerSize', 30, 'Marker', '.', 'LineWidth', 2); %Guardem el objecte plot en una variable per utilitzar les seves propietats mes endavant

range = 1.1 * (l + l); axis([-range range -range range]); axis square;

set(gca, 'nextplot', 'replacechildren'); % Diem que en el seguent plot es pinti a partir d'on acaba l'anterior:

Xcoord = [0, l * sin(phi1), l * sin(phi1) + l * sin(phi2)];
Ycoord = [0, -l * cos(phi1), -l * cos(phi1) - l * cos(phi2)];
set(h, 'XData', Xcoord, 'YData', Ycoord);
drawnow;

succes = 1;