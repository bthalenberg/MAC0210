# Cálculo de erros

Suponha que, usando um algoritmo, computemos uma solução aproximada $\hat{x}$. Se $x$ for
a solução exata, gostaríamos de avaliar o erro absoluto $\|x - \hat{x}\|$, ou o erro relativo

\[
\frac{\| x - \hat{x}\|}{\|x\|}.
\]

Obviamente, se soubessemos esse erro, saberíamos também a resposta exata. Assim, buscamos
um limitante superior para o erro por meio de quantidades que podemos calcular diretamente.
Uma dessas opções é o *residual*:

\[
\hat{r} = b - A\hat{x}
\]

A partir do residual, podemos dizer que

\[

\frac{\| x - \hat{x}\|}{\|x\|}
\leq
\|A\|\|A^{-1}\|  \frac{\|\hat{r}\|}{\|b\|}
\]


Chamamos $\|A\|\|A^{-1}\|$ de número de condição da matriz A, e denotamos $\kappa(A)$.
Ele está relacionado à proximidade da matriz A à singularidade, i.e., o quão próxima
A está de ter um determinante nulo.

Um algoritmo estável é responsável por produzir um residual pequeno. Isso gera erros
aceitavelmente pequenos se a solução do problema é bem condicionada, i.e., se $\kappa(A)$
não for grande demais. Note que pequeno e grande dependem das circunstâncias do problema
(EP do Mascarenhas vs. enviar um foguete pra Lua), mas o conceito qualitativo é geral.

----------------------------------------------------------------------------------------

Uma observação útil: matrizes ortogonais são perfeitamente condicionais. Assim, transformações
ortogonais são especialmente estáveis :-)

----------------------------------------------------------------------------------------
