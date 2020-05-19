function DF = jaco(F, x)
    % This code give you the Jacobian matrix of the function F evaluated in x.
    % The Jacobian matrix is (n,m) meanwhile the sixe of F is n and the size of
    % x is m.

    % F son les n funcions escalars

    % Input: F(x):R^m ---> R^n
    % x: (m x 1)-vector ; F: (n x 1)-vector
    % Output: DF(x) (n x m) Jacobian matrix at x
    % [df1/dx1 ... df1/dxm]
    % [ ...           ... ]
    % [dfn/dx1 ... dfn/dxm]

    % El mateix q que faria la funcio jacobian(F, v) del matlab

    f1 = feval(F, x); m = length(x); n = length(f1);

    h = sqrt(eps); H = eye(m) * h; % eps = epsilon minima maquina.

    DF = zeros(n, m); % Matriu a omplir

    for j = 1:m

        f2 = feval(F, x + H(:, j)); % Seleccionem una fila de H per fer la derivada direccional, fent una variacio epsilon

        DF(:, j) = (f2 - f1) / h; % Apliquem definiicio derivada direccional
    end

end
