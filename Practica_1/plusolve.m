function x = plusolve(L,U,P,b)
n = length(b);

for i=1:1:n
    b([i P(i)]) = b([P(i) i]);
end

y = FS(L, b);
x = BS(U, y);

end