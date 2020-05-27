function df=funBurger(f)
% f ha de ser columna
N = 128;
global D1 D2
df=0.1*D2*f+f.*D1*f;

end