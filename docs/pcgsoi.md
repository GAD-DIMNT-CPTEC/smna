# PCGSOI – Cálculo do Stepsize (`stpcalc.f90`)

> **Status:** Documentação técnica do algoritmo atual implementado no GSI

O cálculo do *stepsize* utilizado na rotina `pcgsoi` pode parecer complexo à primeira vista.
O código foi cuidadosamente escrito para:

* Minimizar erro de arredondamento
* Minimizar recomputação desnecessária
* Maximizar a precisão da estimativa do passo
* Garantir robustez operacional mesmo em regime fortemente não linear

Como consequência, a lógica não é imediatamente transparente.
Este documento descreve exatamente o que a rotina faz.

---

# 1. Visão Geral do Algoritmo

A rotina realiza um **ajuste quadrático iterativo** da função penalidade ao longo da direção de busca atual.

Para cada cálculo de stepsize, o algoritmo:

1. Avalia **pelo menos 3 valores da penalidade** ao longo da direção.
2. Ajusta uma parábola local.
3. Estima o ponto mínimo da parábola.
4. Se necessário, refina a estimativa (até 5 iterações).
5. Aplica lógica de fallback em caso de não convexidade.

---

# 2. Interface Atual da Subrotina

```fortran
subroutine stpcalc(stpinout, sval, sbias, dirx, dval, dbias, diry, &
                   penalty, penaltynew, pjcost, pjcostnew, end_iter)
```

## Entradas

* `stpinout` — estimativa inicial do stepsize
* `sval` — estado atual (representação por bins de observação)
* `sbias` — coeficientes atuais de bias
* `dirx` — direção de busca no espaço de controle
* `dval` — direção de busca nas variáveis de estado
* `dbias` — direção de busca nos parâmetros de bias
* `diry` — aplicação do operador $B^{-1}$ sobre `dirx`
* `penalty` — valor atual da penalidade total
* `pjcost` — componentes principais da penalidade

## Saídas

* `stpinout` — stepsize final estimado
* `penaltynew` — nova penalidade após aplicação do passo
* `pjcostnew` — componentes atualizadas
* `end_iter` — indica falha fatal no cálculo do passo

---

# 3. Formulação Matemática

Ao longo da direção de busca:

$$
x(\alpha) = x_k + \alpha d_k
$$

A penalidade é aproximada localmente por:

$$
J(\alpha) = a + 2b\alpha + c\alpha^2
$$

O mínimo da parábola ocorre em:

$$
\alpha_{min} = -\frac{b}{c}
$$

Para a penalidade total:

$$
\alpha_{min} = -\frac{\sum_i b_i}{\sum_i c_i}
$$

A condição de convexidade local é:

$$
\sum_i c_i > 0
$$

---

# 4. Representação Interna dos Coeficientes

No código, os coeficientes são armazenados como:

```fortran
pstart(1,*) = a
pstart(2,*) = -b
pstart(3,*) = c
```

Portanto, internamente:

$$
J(\alpha) = \text{pstart}(1) - 2,\text{pstart}(2),\alpha + \text{pstart}(3),\alpha^2
$$

Logo:

$$
b = -\text{pstart}(2)
$$

Esse detalhe é importante para interpretar corretamente os sinais.

---

# 5. Termos Tratados como Quadráticos Exatos

Atualmente:

```fortran
npenlin = 3
```

Os três primeiros termos são analiticamente quadráticos ao longo da direção.

## 5.1 Termo de Background

$$
a = \langle x_{save}, y_{save} \rangle
$$

Idealmente:

$$
\langle x_{save}, y_{dir} \rangle =
\langle x_{dir}, y_{save} \rangle
$$

Devido a diferenças de arredondamento, utiliza-se a média.

$$
c = \langle x_{dir}, y_{dir} \rangle
$$

---

## 5.2 Filtro Digital

Calculado em `stpjcdfi`.

---

## 5.3 Restrição de Massa Seca

Calculado em `stpjcpdry`.

---

# 6. Estrutura Diferencial da Matriz `pbc`

A matriz:

```fortran
pbc(i,j)
```

armazenada como:

* `pbc(1,*)` = penalidade em $\alpha_1$
* `pbc(2,*)` = $J(\alpha_2) - J(\alpha_1)$
* `pbc(3,*)` = $J(\alpha_3) - J(\alpha_1)$
* `pbc(4,*)` = $J(0) - J(\alpha_1)$ (primeira iteração)

Ou seja, o código trabalha com:

$$
\Delta J(\alpha) = J(\alpha) - J(\alpha_0)
$$

Isso reduz significativamente cancelamento catastrófico quando:

$$
J(\alpha_2) \approx J(\alpha_1)
$$

---

# 7. Estrutura Detalhada do Array `pbc`

A matriz `pbc` é a estrutura central utilizada para armazenar as contribuições individuais de cada componente da função penalidade ao longo das variações de stepsize avaliadas.

Sua organização é bidimensional:

```fortran
pbc(i,j)
```

onde:

* **Primeira dimensão (`i`)** → variação do stepsize (`sges`)
* **Segunda dimensão (`j`)** → componente específica da penalidade

Essa organização permite:

* Ajuste quadrático independente por componente
* Soma consistente das contribuições globais
* Redução de erro de arredondamento via armazenamento diferencial

---

## 7.1 Primeira Dimensão — Variação do Stepsize

O índice da primeira dimensão representa diferentes valores testados do stepsize:

| Índice     | Significado                                   |
| ---------- | --------------------------------------------- |
| `pbc(1,*)` | Penalidade avaliada em `sges(1)`              |
| `pbc(2,*)` | Penalidade(`sges(2)`) − Penalidade(`sges(1)`) |
| `pbc(3,*)` | Penalidade(`sges(3)`) − Penalidade(`sges(1)`) |
| `pbc(4,*)` | Penalidade(`sges(4)`) − Penalidade(`sges(1)`) |

Ou seja, exceto a primeira linha, os valores armazenados são **diferenças relativas**:

$$
\Delta J(\alpha_k) = J(\alpha_k) - J(\alpha_1)
$$

Essa estratégia reduz cancelamento catastrófico quando:

$$
J(\alpha_k) \approx J(\alpha_1)
$$

e é fundamental para manter precisão em regime de penalidades grandes (ordem 10⁶–10⁹) com reduções pequenas.

---

## 7.2 Segunda Dimensão — Componentes da Penalidade

A segunda dimensão organiza explicitamente as contribuições individuais da função custo:

$$
J = \sum_j J_j
$$

Essa decomposição é essencial para:

* Diagnóstico físico
* Ajuste quadrático componente a componente
* Monitoramento operacional

---

## 7.3 Termos Lineares (Quadráticos Exatos)

Intervalo:

```fortran
ipenlin = 3
pbc(*,1:ipenlin)
```

| Índice     | Componente                                                                |
| ---------- | ------------------------------------------------------------------------- |
| `pbc(*,1)` | Contribuição do **background** + bias de radiância + bias de precipitação |
| `pbc(*,2)` | Termo do **filtro digital**                                               |
| `pbc(*,3)` | Restrição de **pressão seca global** (Jc)                                 |

Esses termos são tratados como **analiticamente quadráticos** ao longo da direção de busca.

Portanto, seus coeficientes $b$ e $c$ são computados exatamente (dentro de erro de máquina).

---

## 7.4 Termos Não Lineares

Intervalo:

```fortran
pbc(*, ipenlin+1 : ipen)
```

Esses termos podem ser fortemente não lineares ao longo da direção.

### Restrições Físicas (Jl / Jc)

| Índice     | Componente                |
| ---------- | ------------------------- |
| `pbc(*,4)` | Umidade negativa (Jl/Jq)  |
| `pbc(*,5)` | Umidade excessiva (Jl/Jq) |
| `pbc(*,6)` | Rajada negativa           |
| `pbc(*,7)` | Visibilidade negativa     |
| `pbc(*,8)` | PBLH negativo             |

Esses termos representam **penalizações físicas de limite inferior ou superior**.

---

### Termos Observacionais (Jo)

Os termos observacionais são organizados explicitamente por tipo:

| Índice      | Observação                 |
| ----------- | -------------------------- |
| `pbc(*,9)`  | Pressão de superfície (ps) |
| `pbc(*,10)` | Temperatura (t)            |
| `pbc(*,11)` | Vento (u/v)                |
| `pbc(*,12)` | Umidade (q)                |
| `pbc(*,13)` | Velocidade do vento (spd)  |
| `pbc(*,14)` | SRW                        |
| `pbc(*,15)` | RW                         |
| `pbc(*,16)` | DW                         |
| `pbc(*,17)` | SST                        |
| `pbc(*,18)` | PW                         |
| `pbc(*,19)` | Precipitação (pcp)         |
| `pbc(*,20)` | Ozônio (oz)                |
| `pbc(*,21)` | O3L (não utilizado)        |
| `pbc(*,22)` | GPS                        |
| `pbc(*,23)` | Radiância                  |
| `pbc(*,24)` | TCP                        |
| `pbc(*,30)` | Rajada observacional       |
| `pbc(*,31)` | Visibilidade observacional |
| `pbc(*,32)` | PBLH observacional         |

⚠ Observação importante:

A posição pode variar conforme novas observações sejam adicionadas ao sistema.
O layout não é conceitualmente fixo — apenas convencional na implementação atual.

---

## 7.5 Interpretação Matemática

Cada coluna `j` da matriz representa:

$$
J_j(\alpha)
$$

A penalidade total é:

$$
J(\alpha) = \sum_j J_j(\alpha)
$$

O ajuste quadrático é feito componente a componente:

$$
J_j(\alpha) \approx a_j + 2b_j \alpha + c_j \alpha^2
$$

O mínimo total é então estimado por:

$$
\alpha_{min} = -\frac{\sum_j b_j}{\sum_j c_j}
$$

A convexidade global exige:

$$
\sum_j c_j > 0
$$

Mesmo que alguns $c_j < 0$ (não linearidade local), a soma total deve ser positiva.

---

## 7.6 Importância Estrutural do `pbc`

O `pbc` permite:

* Separação física das contribuições
* Diagnóstico detalhado de Jb, Jo, Jc, Jl
* Ajuste robusto mesmo com milhões de observações
* Controle de não convexidade local
* Monitoramento da evolução da penalidade por tipo observacional

Sem essa estrutura, seria impossível manter:

* Estabilidade numérica
* Transparência diagnóstica
* Robustez operacional

---

# 8. Ajuste Quadrático

Para cada componente são estimados:

$$
b_i, c_i
$$

O mínimo global estimado é:

$$
\alpha_{min} = -\frac{\sum_i b_i}{\sum_i c_i}
$$

Pode haver termos individuais com $c_i < 0$,
mas exige-se:

$$
\sum_i c_i > 10^{-20}
$$

---

# 9. Processo Iterativo

Máximo de iterações:

```fortran
istp_iter = 5
```

Em cada iteração:

```fortran
sges(1) = melhor anterior
sges(2) = (1 - dels) * melhor
sges(3) = (1 + dels) * melhor
```

`dels` inicia em 0.1 e pode diminuir.

---

# 10. Critério de Convergência

$$
\left|\frac{\Delta J}{J}\right| < 10^{-17}
$$

Indica convergência extremamente restritiva.

---

# 11. Casos Patológicos

## 11.1 Quadrática Não Convexa

Se:

$$
\sum_i c_i \le 10^{-20}
$$

ou passo negativo:

1. Verifica se algum passo positivo reduz penalidade.
2. Seleciona o de maior redução.
3. Caso contrário:

$$
\alpha = 0.1 \times \min(\alpha_{positivo})
$$

---

# 12. Condição Fatal

Se o passo final permanece negativo:

* Mensagens diagnósticas são impressas.
* `end_iter = .true.`
* Minimização interna é interrompida.

---

# 13. Atualização Final

Se passo válido encontrado:

* Atualiza estado.
* Atualiza bias.
* Recalcula `penaltynew`.
* Atualiza `pjcostnew`.

---

# 14. Filosofia Numérica

A rotina assume que:

* O termo de background é estritamente quadrático.
* O termo observacional pode ser fortemente não linear.
* A paisagem da penalidade pode não ser perfeitamente convexa.

Para garantir robustez operacional:

* Usa armazenamento diferencial.
* Usa curvatura analítica para termos lineares.
* Usa diferenças finitas para termos não lineares.
* Implementa fallback protegido.
* Impõe limite mínimo de curvatura.

Isso permite estabilidade mesmo com milhões de observações.

---

# 15. Resumo Final

A rotina `stpcalc`:

* Ajusta uma quadrática local à penalidade total.
* Estima o stepsize ótimo ao longo da direção de busca.
* Trata termos lineares analiticamente.
* Trata termos não lineares via diferenças finitas.
* Minimiza erro de arredondamento.
* Aplica salvaguardas robustas contra não convexidade.

É um componente central da robustez do PCGSOI no GSI.

