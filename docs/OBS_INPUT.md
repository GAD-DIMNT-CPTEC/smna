# GSI – Análise Técnica Definitiva de `dval`, `dthin` e `dsfcalc`
## Sistema: Gridpoint Statistical Interpolation (GSI)
## Aplicação: SMNA
## Base: Inspeção direta do código-fonte (NCEP/GSI)

---

# 1. Objetivo

Este documento estabelece a interpretação técnica correta dos parâmetros:

- `dval`
- `dthin`
- `dsfcalc`

no bloco `OBS_INPUT::` do GSI, com base em análise direta do código-fonte e verificação contra a configuração operacional do NCEP.

O objetivo é eliminar interpretações simplificadas e documentar o comportamento real no sistema.

---

# 2. Estrutura do `OBS_INPUT::`

Cada linha define um grupo de observações:

```

OBS_INPUT::
dfile   dtype   dplat   dsis   dval   dthin   dsfcalc

````

O parsing ocorre em:

- `obsmod.F90 → init_instr_table_`

Armazenamento interno:

```fortran
real(r_kind), allocatable :: dval(:)
integer, allocatable :: dthin(:)
integer, allocatable :: dsfcalc(:)
````

Um flag global é definido:

```fortran
if (dval(ii) > 0.0) dval_use = .true.
```

Este detalhe é fundamental para a interpretação correta.

---

# 3. `dval` – Peso Relativo Opcional no Thinning

## 3.1 Papel real

`dval` é um **fator opcional de ponderação relativa dentro do algoritmo de thinning/superob para dados de satelite**.

Ele não é obrigatório.
Ele não é usado no modo padrão.
Ele só entra em ação se algum valor for maior que zero.

---

## 3.2 Modo padrão (todos `dval = 0.0`)

Se todas as entradas no `OBS_INPUT` tiverem:

```
dval = 0.0
```

Então:

* `dval_use = .false.`
* O algoritmo de thinning ignora completamente o weighting por `dval`
* O sistema opera em modo thinning padrão

Isso é exatamente o que ocorre no NCEP operacional.

Portanto:

> `dval = 0.0` NÃO desativa observação.
> `dval = 0.0` NÃO impede assimilação.
> `dval = 0.0` é o comportamento padrão.

---

## 3.3 Modo ponderado (algum `dval > 0.0`)

Se ao menos um grupo tiver `dval > 0.0`:

* `dval_use = .true.`
* O algoritmo de `satthin` passa a calcular:

$$
ratio_i = \frac{dval_i}{\sum_j dval_j}
$$

Esta razão controla qual observação sobrevive dentro do grid box.

Esse é um modo especial de operação raramente usado no NCEP.

---

## 3.4 Conclusão técnica sobre `dval`

| Situação           | Comportamento                           |
| ------------------ | --------------------------------------- |
| Todos `dval = 0.0` | Thinning padrão (modo operacional NCEP) |
| Algum `dval > 0.0` | Ativa thinning ponderado                |
| Valores diferentes | Competição relativa intra-grid          |

---

# 4. `dthin` – Controle do Thinning

## 4.1 Papel

`dthin` controla se o algoritmo de thinning será aplicado.

Ele é passado como:

```fortran
ithin = dthin(i)
```

e atua diretamente em `satthin.F90`.

---

## 4.2 Interpretação prática

| Valor | Comportamento                                  |
| ----- | ---------------------------------------------- |
| -1    | Sem thinning                                   |
| 0     | Thinning padrão                                |
| 1     | Thinning ativado                               |
| 2,3   | Modos especiais (dependentes da versão/sensor) |

Se `dthin = -1`, o thinning não ocorre e `dval` perde qualquer relevância.

---

# 5. `dsfcalc` – Recalcular Pressão de Superfície

## 5.1 Papel

`dsfcalc` controla se o GSI recalcula a pressão de superfície (`psfc`) para radiâncias.

Ele é passado aos readers de satelites e influencia o operador radiativo (CRTM).

---

## 5.2 Interpretação prática

| Valor | Comportamento             |
| ----- | ------------------------- |
| 0     | Não recalcula psfc        |
| 1     | Recalcula psfc com modelo |

Impacta especialmente sensores sensíveis à estrutura vertical.

---

# 6. Síntese Conceitual

| Parâmetro | Atua onde    | Papel Real                      |
| --------- | ------------ | ------------------------------- |
| `dval`    | satthin      | Ponderação opcional no thinning |
| `dthin`   | satthin      | Ativa/desativa thinning         |
| `dsfcalc` | readers/CRTM | Recalcula psfc                  |

Todos atuam antes da minimização variacional.

Nenhum altera diretamente R ou a função custo.

---

# 7. Decisão Operacional – Manual Interno SMNA

## Caso A – Configuração padrão operacional

```
dval = 0.0
dthin >= 0
```

→ Thinning padrão
→ Comportamento idêntico ao NCEP

---

## Caso B – Ativar ponderação relativa entre sensores

```
dval > 0.0 (valores diferentes)
dthin >= 0
```

→ Ativa modo ponderado
→ Sensor com maior `dval` tem maior chance de sobreviver

---

## Caso C – Sem thinning

```
dthin = -1
```

→ Todas as observações entram
→ `dval` torna-se irrelevante

---

## Caso D – Teste físico vertical

```
dsfcalc = 0 vs 1
```

→ Avaliar impacto em sensores micro-ondas

---

# 8. Observação Crítica

A tabela operacional do NCEP demonstra claramente:

```
dval = 0.0 para todos os sensores
```

Logo:

> O sistema operacional padrão NÃO utiliza weighting por `dval`.

Esse parâmetro é um recurso avançado, não um controle básico.

---

# 9. Conclusão Final

* `dval` é um mecanismo opcional de ponderação no thinning.
* O padrão operacional é mantê-lo zerado.
* Ele não desativa observações.
* Ele só altera comportamento se explicitamente ativado.
* `dthin` é o verdadeiro controlador do thinning.
* `dsfcalc` atua no ajuste vertical do operador radiativo.

---

Documento validado contra inspeção direta de:

* `obsmod.F90`
* `read_obs.F90`
* `satthin.F90`
* `read_bufrtovs.f90`
* `read_atms.f90`
* `read_satwnd.f90`
* `read_gps.f90`
* `crtm_interface.f90`

---

