# Fatoração QR

A fatoração QR decompõe uma matriz A em uma matriz Q com colunas ortonormais e uma
matriz R triangular superior. A deve ter colunas linearmente independentes. Caso
A seja quadrada, então $Q^{-1} = Q^T$.

Para quaisquer matriz ortogonal $P$ e vetor $w$ podemos escrever

\[
\|b - Aw\| = \|Pb - PAw\|.
\]

Nosso objetivo é obter uma transformação tal que

\[
A = Q \begin{pmatrix} R \\ 0 \end{pmatrix}
\]

## Motivação

Esse tipo de fatoração é especialmente útil para resolver o problema dos mínimos
quadrados (não me pergunte como e nem por quê), assim como o algoritmo de SVD. O método das equações normais, da forma $B = A^TA$, que usa a fatoração de Cholesky, pode ter $\kappa(A)^2$. A intenção do algoritmo QR, então, é minimizar esse erro,
obtendo $\kappa(A)$. Ele realiza mais operações -- se implementado de maneira eficiente, faz em média $2mn^2 - \frac{2}{3}n^3$ flops -- mas é mais robusto.

## A transformação Householder

A maneira mais fácil de realizar isso é a ortogonalização de Gram Schmidt, que já
vimos em Álgebra Linear. Obviamente, como nada na vida é simples, não é esse o método
que precisamos saber para a prova. Vamos usar a transformação Householder, porque
aparentemente ela é melhor para propósitos gerais e a mais semelhante à fatoração LU.

Suponha que dado um vetor z, quiséssemos encontrar uma transformação ortogonal que
zerasse todas as entradas do vetor exceto pela primeira. Uma transformação possível é
$Pz = z - 2uu^Tz$, onde u é um vetor unitário e I é uma matriz identidade. P, na sua
forma geal, $P = I - 2uu^T$ é chamado de refletor, uma vez que $Pu = -u$. No caso,
o vetor unitário será na direção de $z \pm \|z\|e_1$. O sinal será selecionado de acordo
com $z_1$, para reduzir a possibilidade de erros de cancelamento.

Esse processo pode ser aplicado e extendido de maneira análoga à derivação da fatoração LU.

## O algoritmo

Para realizar a fatoração QR, procedemos da seguinte maneira:

1. Aplique a reflexão $P^{(1)}$ que transforma a primeira coluna de A em um múltiplo de $e_1$, usando a primeira coluna de A no lugar de z descrito acima.
2. Em seguida, aplique uma transformação semelhante na segunda coluna de $P^{(1)}A$,
zerando os elementos abaixo da diagonal.
3. Continue realizando transformações análogas até obter uma matriz triangular superior.

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

Código em Java (da biblioteca Jama):
```
public QRDecomposition (Matrix A) {
      // Initialize.
      QR = A.getArrayCopy();
      m = A.getRowDimension();
      n = A.getColumnDimension();
      Rdiag = new double[n];

      // Main loop.
      for (int k = 0; k < n; k++) {
         // Compute 2-norm of k-th column without under/overflow.
         double nrm = 0;
         for (int i = k; i < m; i++) {
            nrm = Maths.hypot(nrm,QR[i][k]);
         }

         if (nrm != 0.0) {
            // Form k-th Householder vector.
            if (QR[k][k] < 0) {
               nrm = -nrm;
            }
            for (int i = k; i < m; i++) {
               QR[i][k] /= nrm;
            }
            QR[k][k] += 1.0;

            // Apply transformation to remaining columns.
            for (int j = k+1; j < n; j++) {
               double s = 0.0; 
               for (int i = k; i < m; i++) {
                  s += QR[i][k]*QR[i][j];
               }
               s = -s/QR[k][k];
               for (int i = k; i < m; i++) {
                  QR[i][j] += s*QR[i][k];
               }
            }
         }
         Rdiag[k] = -nrm;
      }
   }
```