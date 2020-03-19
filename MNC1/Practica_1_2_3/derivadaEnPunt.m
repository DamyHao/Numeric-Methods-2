function out = derivadaEnPunt(F, x)
%Amb els més petits que menys 16 no funciona (surt 0).
%Amb alguns nombres entre -10 y -16 te un comportament extrany ja que
%sembla tenir menys precisió.
minNumber = 10^-10;
out = (F(x + minNumber) - F(x)) / minNumber;

end