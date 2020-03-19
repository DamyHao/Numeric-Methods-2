function [Ximp,it] = newton(ini,fin,tol,itmax,fun) %b
it=0; Ximp=[(ini+fin),(ini+fin)/2,(ini+fin)/3]; tolk=1;
  if fun(ini)*fun(fin)<0 
      while tolk>tol && it<itmax
          Ximp = [Ximp Ximp(end)-(fun(Ximp(end))/der(Ximp(end),fun))];
          it = it+1;
          tolk = abs(Ximp(end)-Ximp(end-1));
      end
  else
      disp('Les dues imatges tenen el mateix signe.');
  end
end

%[xk,it]=newton(60,120,10^-10,15,@diferencia_newton)
%considero ini=60 i fin=120 donat que el punt de tall que m'interessa se situa entre
%aquests dos valors, tal com es pot veure a simple vista a la gràfica de
%l'apartat anterior.