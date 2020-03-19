function superficie = superficie_cilindre (radi, altura);
area_base = pi*(radi^2);
circumferencia = 2*pi*radi;
superficie = 2*(area_base)*circumferencia*altura;

%clear all;
%radi=2; altura=10;
%area=superficie_cilindre(radi, altura);
%area