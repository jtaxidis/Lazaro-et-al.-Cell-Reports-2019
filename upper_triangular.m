function A = upper_triangular(A)

if all(size(A) > 1)
    A = A + 1000;
    A = triu(A,1);
    A = nonzeros(A);
    A = A - 1000;
end