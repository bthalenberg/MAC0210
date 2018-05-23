# Decomposição LU

O método de eliminação de Gauss Jordan decompõe uma matriz A em um produto LU de uma
matriz triangular inferior L e uma matriz triangular superior U. Uma matriz é triangular
superior se $a_{ij} = 0 \forall i > j$ e triangular inferior se $a_{ij} = 0 \forall i < j$.
Especificamente, a matriz L guarda os multiplicadores $l_{ik} = \frac{a_{ik}}{a_{kk}}$,
enquanto a matriz U guarda o resultado das operações sobre A, ou seja, a matriz obtida ao
final da fase de substituição.

Para resolver uma equação Ax = b por decomposição LU, fazemos:

```
Obtenha A = LU
Resolva Ly = b (substituição)
Resolva Ux = y (retrossubstituição)
```

## Motivação

Tome duas equações com a mesma matriz, mas com lados direitos diferentes:

\[
Ax = b
\]
\[
Ax' = b'
\]

Caso fossemos resolver as duas isoladamente, teríamos um custo de $2(f + s)$, onde
$f$ é o número de operações para fatorar A e $s$ o número de operações para resolver
o sistema. Ao fatorar A de maneira aproveitável para ambas, temos um custo de $f + 2s$.
Se $f \gg s$, resolvemos duas equações pelo custo de uma.

A decomposição LU também encontra aplicações na resolução de determinantes, pois o
determinante de um produto de matrizes é igual ao produto do determinante das matrizes,
e o determinante de uma matriz triangular é o produto dos elementos da diagonal principal. Assim, $det(A) = det(L).det(U) = \prod_{i=0}{n} u_ii$.

Podemos, ainda, utilizar a fatoração para inverter matrizes. Definimos $e_i$ como
os vetores coordenada padrões de $\Re^n$. Para cada um deles, resolvemos
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

### O algoritmo

```
void lu(float *m, int n, int *p) {
    for (int i = 0; i < n; i++) p[i] = i;
    for (int i = 0; o < n - 1; i++) {
        int i_max = 0;
        float v_max = fabs(*m);
        float *c = m;
        // pivotamento parcial
        // procuramos o maior item da coluna
        for (int k = i + 1; i < n; i++) {
            if (fabs(*(c += n)) > v_max) {
                i_max = i;
                v_max = fabs(*c);
            }
        }
        // fazemos a troca de pivô
        if (i_max > i) {
            float *pa = m - i;
            float *pb = pa + (i_max - i) % n;
            for (int j = 0; j < n; j++; pa++; pb++)
                swap(*pa, *pb);
            swap(p[i], p[i_max]);
        }
        float **pc = m;
        // fatoração
        for (int j = i + 1; j < n; j++) {
            **(pc += n) /= *m;
            float *pa = pc;
            float *pb = m;
            for (int k = i + 1; k < n; k++)
                **(++pa) -= *pc * *(++pa);
                m += n + 1;
        }
    }
}
```
