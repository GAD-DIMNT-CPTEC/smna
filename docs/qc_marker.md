# 📘 Marcas de Controle de Qualidade (QC Markers) no PREPBUFR e GSI

Esta documentação descreve o funcionamento das **marcas de controle de qualidade (QC markers)** utilizadas no sistema **PREPBUFR → GSI**, com base:

* No comportamento observado no código Fortran
* Na lógica implementada em `statcount`
* No fluxo operacional do GSI

---

# 1️⃣ Conceito Geral

As **QC markers** são valores inteiros associados a cada observação no arquivo PREPBUFR.

Elas indicam:

* A qualidade da observação
* Se passou nas verificações preliminares
* Se está apta para assimilação
* Se será apenas monitorada
* Se deve ser descartada

O GSI utiliza essas marcas para decidir o destino final da observação.

---

# 2️⃣ Variáveis Fundamentais no Processo de QC

Para interpretar corretamente o controle de qualidade, devemos sempre analisar:

| Variável  | Papel                            |
| --------- | -------------------------------- |
| `idqc`    | Marca de qualidade (PREPBUFR)    |
| `iuse`    | Indicador final de uso no GSI    |
| `noiqc`   | Parâmetro do namelist            |
| `lim_qm`  | Limite entre monitoramento e uso |
| `varName` | Variável (ps, uv, t, q...)       |
| `varType` | Tipo de instrumento (KX)         |

---

# 3️⃣ Interpretação da Marca `idqc`

A tabela operacional é:

| Intervalo de `idqc`             | Processo no GSI                  |
| ------------------------------- | -------------------------------- |
| `idqc > 15` ou `idqc <= 0`      | Observação descartada na leitura |
| `idqc >= lim_qm` e `idqc <= 15` | Observação monitorada            |
| `idqc > 0` e `idqc < lim_qm`    | Elegível para assimilação        |

---

# 4️⃣ Papel do parâmetro `noiqc`

O valor de `lim_qm` depende da opção `noiqc` no namelist do GSI.

## 🔹 Quando `noiqc = True` (sem OI QC)

| Variável | lim_qm |
| -------- | ------ |
| ps       | 7      |
| outras   | 8      |

---

## 🔹 Quando `noiqc = False` (com OI QC)

| Variável | lim_qm |
| -------- | ------ |
| todas    | 4      |

---

# 5️⃣ Papel da variável `iuse`

`iuse` representa o estado final após todas as verificações internas do GSI.

| Valor | Significado                 |
| ----- | --------------------------- |
| `1`   | Observação usada na análise |
| `-1`  | Observação não usada        |

⚠ Importante:
Uma observação pode ter `idqc` dentro da faixa aceitável e ainda assim ser rejeitada posteriormente.

---

# 6️⃣ Classificação Final

A lógica implementada no código é:

### 🔹 Assimiladas

```
iuse == 1
```

---

### 🔹 Monitoradas

```
(iuse == -1) AND (idqc >= lim_qm AND idqc <= 15)
```

---

### 🔹 Rejeitadas

```
(iuse == -1) AND (
    idqc > 15 OR
    idqc <= 0 OR
    idqc < lim_qm
)
```

---

# 7️⃣ Fluxo Completo do Processo de QC

```
Observação
    ↓
QC no PREPBUFR
    ↓
Atribuição de idqc
    ↓
Leitura no GSI
    ↓
Aplicação de lim_qm (dependente de noiqc)
    ↓
Gross check / inovação / background check
    ↓
Definição de iuse
    ↓
Classificação final:
    - Assimilada
    - Monitorada
    - Rejeitada
```

---

# 8️⃣ O que cada classe significa fisicamente

| Classe     | Interpretação         |
| ---------- | --------------------- |
| Assimilada | Contribuiu na análise |
| Monitorada | Lida, mas não usada   |
| Rejeitada  | Descartada            |

---

# 9️⃣ Relação com o arquivo `convinfo`

O `convinfo` define:

* Erro observacional (σ)
* Limiar do gross check
* Parâmetros de thinning
* Parâmetros VarQC

Esses parâmetros influenciam:

* Se a observação será rejeitada
* Se será usada
* O valor final de `iuse`

---

# 🔟 Relação com arquivos diag

Nos arquivos diag:

* `idqc` representa o estado herdado do PREPBUFR
* `iuse` representa o estado final após o GSI

Para análise estatística correta, deve-se sempre usar **ambos**.

---

# 1️⃣1️⃣ Como verificar QC corretamente

## Para saber se foi usada

```python
iuse == 1
```

---

## Para saber se foi monitorada

```python
(iuse == -1) & (idqc >= lim_qm)
```

---

## Para saber se foi rejeitada

```python
iuse == -1
```

mas separando rejeição física de monitoramento.

---

# 1️⃣2️⃣ Resumo Conceitual

| Etapa                  | Variável relevante |
| ---------------------- | ------------------ |
| QC inicial             | idqc               |
| Elegibilidade          | idqc + lim_qm      |
| Uso final              | iuse               |
| Estatística científica | idqc + iuse        |

---

# 1️⃣3️⃣ Observação Científica Importante

* `idqc` representa qualidade prévia
* `iuse` representa decisão final
* `noiqc` altera o limiar
* `lim_qm` muda por variável
* O GSI pode rejeitar dados mesmo com `idqc` válido

---

# 1️⃣4️⃣ Conclusão

As marcas de controle de qualidade no sistema PREPBUFR + GSI formam um sistema hierárquico composto por:

* Avaliação prévia (idqc)
* Aplicação de limiar (lim_qm)
* Verificações internas do GSI
* Decisão final de uso (iuse)

A análise correta do comportamento da assimilação exige considerar simultaneamente:

```
idqc
iuse
noiqc
lim_qm
varName
varType
```