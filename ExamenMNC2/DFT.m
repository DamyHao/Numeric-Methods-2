function fk = DFT(fj)
    % Li donem els punts avaluats en physical space i ens retorna els
    % coeficients de la serie de Fourier (transformades de fourier (fourier space))
    % tot en columna
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
    
    N = length(fj);
    W = exp(-1i * 2 * pi / N);
    k = 0:N - 1; j = 0:N - 1;

    F = (1 / N) * W.^((k' - N / 2) * j);
    fk = F * fj;

end
