function B = baricentrica(z, xj,fun)
% funz funcio avaluada en z

numerador = [];
for j=0:1:length(xj)+1
    numerador = numerador + (-1)^j.*(z-xj(j+1))^-1;
    
    denominador = (-1)^j.*(z-xj(j+1))^-1;
    
end;

end