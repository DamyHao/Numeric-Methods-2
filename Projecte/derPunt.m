function valor = derPunt(x,fun)
%   Primera derivada en un cert punt x de la funcio fun
    H = 10^-15;
    valor = (fun(x)+fun(x+h))/H;
end

