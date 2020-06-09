function unique = uniqueSol(solucions)
% Punts son columnes
epsilon = sqrt(eps);
[~, n] = size(solucions);
unique = [solucions(:,1)];
nsol = 1;
for ii=1:1:n
    jj = 1;
    goOut = 1;
    while jj <= nsol && goOut
        %Solucio repetida
        if( max(abs(solucions(:, ii) - unique(:, jj))) < epsilon)
            goOut = 0;
        end
        jj = jj + 1;
    end
    if(goOut == 1)
        unique = [unique, solucions(:, ii)];
        nsol = nsol + 1;
    end
end
end