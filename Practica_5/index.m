%Provant...
%disp(genpath([fileparts(pwd), "\Practica_1", 'PLU.m']))
%addpath(genpath(fileparts(pwd), "\Practica_2"))
%load(fullfile('..', 'Practica_1', 'PLU.m'))
addpath("../Practica_1")





f = @(x)([exp(-20 * (x(1) + 1/4).^2) + (sin(30 * x(1)) .* exp(-20 * (x(1) - 1/4).^2)) / 4; 1]);

disp(f([1,2]))