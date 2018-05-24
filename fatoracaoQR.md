# Fatoração QR


## O algoritmo

Em MatLab/Octave porque não sou obrigada e muito menos paga pra isso:

```
[m,n] = size(A);
Q = eye(m);
% lembrando que em MatLab índices começam em 1
for k = 1:n
    // Find the HH reflector
    z = A(k:m, k);
    v = [- sign(z(1)) * norm(z) - z(1); z(2:end)];
    v = v / sqrt(v' * v);
    // Apply the HH reflection to each column of A and Q
    for j = 1:n
        A(k:m, j) -= v * 2* (v'A(k:m, j));
    end
    for j = 1:m
        Q(k:m, j) -= v * 2* (v'Q(k:m, j));
    end
 end
 Q = Q';
 R = triu(A); // upper triangular part of A
```

A fatoração QR é mais estável do que a fatoração de Cholesky, o que é importante
quando as colunas são 'quase' linearmente dependentes. Por outro lado, ela é
duas vezes mais lenta.
