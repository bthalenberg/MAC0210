# Decomposição de Cholesky

Uma matriz $A$ é simétrica se $A^T = A$ e definida positiva de $x^TAx > 0\:
\forall x \neq 0$. A decomposição de tais matrizes é relativamente mais simples do que
a decomposição de matrizes comuns, de forma que podemos simplificar o algoritmo da
decomposição LU. Por exemplo, temos que não é necessário fazer pivotamento para
estabilizar o algoritmo.

Considere a decomposição de uma de tais matrizes como $A = LU$. Podemos escrever $U$
como $DU'$, onde D é sua diagonal, e U' tem os elementos da diagonal substituídos por
1s.

Temos, portanto, que $A = LDU'$. Transpondo, obtemos $A = A^T = U'^{T}DL^{T}. No entanto,
como a decomposição é única, devemos ter $U'^{T} = L$. Logo,

\[
A = LDL^T
\]

onde L é triangular inferior unitária, ou seja, com elementos da diagonal iguais a 1,
 e D é diagonal com elementos $u_{kk}$. Também podemos escrever $D = D^{\frac{1}{2}}D^{\frac{1}{2}}$,
com $D^{\frac{1}{2}} = diag{\sqrt{u_{11},...,\sqrt{u_{nn}}}}$, de forma que

\[
A = GG^T
\]

com $G = LD^{\frac{1}{2}}$ triangular inferior.

## Motivação

O algoritmo para a decomposição de Cholesky necessita cerca de metade do espaço requerido
pela decomposição LU, por conta da simetria, e cerca de metade dos flops de LU,
$\frac{1}{3}n^3 + O(n^2)$. Isso pode ser feito porque os cálculos para $A = GG^T$
podem ser feitos elemento a elemento. Note que isso somente ocorre porque os argumentos
da raiz quadrada sempre serão positivos no caso de A simétrica positiva definida,
mas não para matrizes genéricas.

## Algoritmo

```
// Note que as entradas de A são sobrescritas por G
for k = 1:n - 1
    a_kk = sqrt(a_kk)
    for i = k + 1:n
        a_ik = a_ik/a_kk
    end
    for j = k + 1:n
        for i = j:n
            a_ij = a_ij - a_ik * a_jk
        end
    end
end
a_nn = sqrt(a_nn)
```
