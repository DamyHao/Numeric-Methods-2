function nousPunts = nextPoint(anticsPunts, Fz, Fy)
%   Trobar el següent punt equipotencial (anticsPunts, Fz, Fy)
epsi = 0.001;
arrel = sqrt(Fy^2+Fz^2);
% A LA FITXA ESTAN AL REVÉS!!! ;(
f = [-Fy/arrel, Fz/arrel];
nousPunts = anticsPunts + epsi.*f;
end

