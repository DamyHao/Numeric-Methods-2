%% CANVIS QUE HAS DE FER:

AL VOLTANT DE LA LINIA 152:

%Y proporciono los valores que voy obteniendo para cada valor de theta:
display = sprintf('Para theta = %d, el valor de K, según el método utilizado, es:', theta);
disp(display);
I_trap_c
I_cc
I_fej


%%

A LA PRIMERA GRÀFICA:

title('Visualización de la función (error relativo(/Theta) - 0.05)')

%%

A LA SEGONA GRÀFICA:

xlabel('Ángulo (/Theta_0)');
ylabel('K(/Theta_0)');

%%

A LA TERCERA GRÀFICA:

xlabel('Ángulo (/Theta_0)');
ylabel('abs(K_/Theta_f_e_j - K_/Theta_c_c');
