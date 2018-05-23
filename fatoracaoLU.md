# Métodos numéricos de álgebra linear

## Decomposição LU

O método de eliminação de Gauss Jordan decompõe uma matriz A em um produto LU de uma
matriz triangular inferior L e uma matriz triangular superior U. Uma matriz é triangular
superior se $a_{ij} = 0 \forall i > j$ e triangular inferior se $a_{ij} = 0 \forall i < j$.
Especificamente, a matriz L guarda os multiplicadores $ l_{ik} = \frac{a_{ik}}{a_{kk}},
enquanto a matriz U guarda o resultado das operações sobre A, ou seja, a matriz obtida ao
final da fase de substituição.

Para resolver uma equação Ax = b por decomposição LU, fazemos:

```
Obtenha A = LU
Resolva Ly = b (substituição)
Resolva Ux = y (retrossubstituição)
```

### Motivação

Tome duas equações com a mesma matriz, mas com lados direitos diferentes:

$
Ax = b
Ax' = b'
$

Caso fossemos resolver as duas isoladamente, teríamos um custo de $2(f + s)$, onde
$f$ é o número de operações para fatorar A e $s$ o número de operações para resolver
o sistema. Ao fatorar A de maneira aproveitável para ambas, temos um custo de $f + 2s$.
Se $f \gg s$, resolvemos duas equações pelo custo de uma.

A decomposição LU também encontra aplicações na resolução de determinantes, pois o
determinante de um produto de matrizes é igual ao produto do determinante das matrizes,
e o determinante de uma matriz triangular é o produto dos elementos da diagonal principal. Assim, $det(A) = det(L).det(U) = \prod_{i=0}{n} u_ii$.

Podemos, ainda, utilizar a fatoração para inverter matrizes. Definimos $e_i$ como
os vetores coordenada padrões de $\mathds{R}^n$. Para cada um deles, resolvemos
$Ax_i = e_i$, ou, equivalentemente, $L(U x_i) = e_i$. O vetor $x_i$ obtido será
a i-ésima coluna de $A^{-1}$. A inversão toma $\frac{8}{3}n^3$ operações.

Note que é possível armazenar L e U na mesma matriz, uma vez que as suas partes
não-nulas são complementares.

## Decomposição PA = LU

A decomposição LU não pode ser aplicada para todas as matrizes não singulares,
pois não lida bem com pivôs não-nulos. Além disso, erros indesejáveis de
arredondamento acumulado podem surgir a partir de pivôs pequenos.

A fatoração PA = LU visa solucionar esse problema por meio do pivotamento
parcial. Antes de realizar as operações da fase de substituição de uma linha,
busca-se o maior elemento da coluna sendo examinada, e invertem-se as linhas
para que ele torne-se o pivô. A operação de troca de linhas será registrada em
uma matriz de permutações P. Começamos com P = I e, a cada permutação em A,
permutamos as linhas equivalentes de P.

Ao total, são realizadas $\frac{2}{3}n^3$ flops (_floating-point operations_),
mais $n^2$ em cada fase de substituição para resolução do sistema.

```
Inicialize a matriz de permutações P
Para cada linha k:
    Encontre a linha q cujo pivô é o máximo relativo da coluna k
    Permute a linha k com a linha q
    Compute a coluna correspondente de L (l_jk = a_jk/akk)
    Atualize a submatriz

Para cada linha k, 1 <= k < n:
    y(k) = b(k) - A
```

Note que isso ainda não garante a estabilidade para erros de arredondamento.
Para isso, deve ser feito o escalamento do pivotamento parcial. No geral,
o custo computacional é muito grande, então não fazemos isso.

## Fatoração QR

[m,n] = size(A);
Q = eye(m);
for (int k = 0; k < n; k++)
    // Find the HH reflector
    z = A(k:m, k);
    v = [- sign(z(0)) * norm(z) - z(0); z(1:end)];
    v = v / sqrt(v' * v);
    // Apply the HH reflection to each column of A and Q
    for (int j = 0; j < n; j++)
        A(k:m, j) -= v * 2* (v'A(k:m, j));
    for (int j = 0; j < m; j++)
        Q(k:m, j) -= v * 2* (v'Q(k:m, j));
 Q = Q';
 R = triu(A); // upper triangular part of A
