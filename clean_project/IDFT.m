function fj= IDFT (fk)
    % Et dona la funcio en el physicial space evaluada en els punts xj a partir
    % dels coeficients de fourier (en format columna)
    % Aquest es en Ricard, un mapache que te pànic a les N imparells
    %                  __        .-.
    %              .-"` .`'.    /\\|
    %      _(\-/)_" ,  .   ,\  /\\\/
    %     {(#b^d#)} .   ./,  |/\\\/
    %     `-.(Y).-`  ,  |  , |\.-`
    %          /~/,_/~~~\,__.-`
    %         ////~    // ~\\
    %       ==`==`   ==`   ==`
    % No espantis al Ricard, posa N parells

N = length(fk); WN = exp(-i*2*pi/N); jj = [0:N-1]'; kk = jj'; 
P = WN.^(-jj*(kk-N/2)); fj = P*fk

end