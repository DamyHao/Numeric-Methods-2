clear all
dni1 = [3 9 4 2 3 0 1 5]
dni2 = [4 7 5 9 3 8 4 5]
z = [2 5 7]
dni1(z);
dni2(z);
dni1_medio = sum(dni1)/length(dni1);
dni2_medio = sum(dni2)/length(dni2);
suma_medio = sum(dni1+dni2)/length(dni1+dni2);
mul = dni1.*dni2;
mul_orden = sort(mul);
mul_orden_trozo = [mul_orden(1:3) mul_orden(6:end)];
dni1_orden = fliplr(sort(dni1));
dni2_orden = fliplr(sort(dni2));
dni1_orde_trozo = dni1_orden(z);
dni2_orden_trozo