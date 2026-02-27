# Ajuste do DeltaT no BAM/SMNA

## 🧩 Contexto

Durante os testes de estabilidade do **SMNA** (Sistema de Modelagem Numérica e Assimilação) com o modelo **BAM** em truncamento **T299L64**, foi identificada uma **falha de segmentação (`segmentation fault`)** após aproximadamente **15.000 segundos de simulação**. A investigação indicou que a origem do problema estava relacionada à configuração do **timestep (`DeltaT`)**.

---

## 🔍 Diagnóstico

A análise inicial apontou que a instabilidade poderia ser causada por um ou mais dos seguintes fatores:

* **Timestep próximo do limite de estabilidade**, amplificando oscilações numéricas em modos de alta frequência.
* **Acúmulo de erros numéricos** nos arquivos de *restart*, não totalmente corrigidos ao longo dos ciclos de assimilação.
* **Violação da condição CFL (Courant–Friedrichs–Lewy)**, associada a variações dinâmicas intensas na atmosfera.

Como primeira ação, foi realizada a **reinicialização do modelo com `ColdStart` (initlz = 2)**, eliminando possíveis valores espúrios acumulados. No entanto, o problema persistiu.

---

## ⚙️ Soluções Testadas

| Teste | Configuração                | Resultado                                           |
| ----- | --------------------------- | --------------------------------------------------- |
| **1** | `DeltaT = 240 s` (original) | Instabilidade e *segmentation fault* após ~15.000 s |
| **2** | `DeltaT = 220 s`            | Integração concluída sem falhas (*estável*)         |
| **3** | `DeltaT = 225 s`            | Persistência do erro                                |
| **4** | `DeltaT = 200 s`            | Integração estável e completa                       |

Esses resultados indicam que o problema não estava exclusivamente ligado à relação de submúltiplos de 3600 s, mas também à sensibilidade do modelo sob determinadas condições dinâmicas.

---

## 🧠 Considerações de Bonatti

> O `DeltaT` deve ser **submúltiplo de 3600 s**, garantindo chamadas corretas à radiação em horas inteiras.
>
> Nas versões anteriores:
>
> * **Radiação de onda longa:** chamada a cada **3 horas**
> * **Radiação de onda curta:** chamada a cada **1 hora**
>
> Exemplo de submúltiplos válidos:
> `225 s (3600/16)` • `240 s (3600/15)` • `200 s (3600/18)`

Além disso, o BAM contém um **sistema interno de controle do CFL** (filtragem de ondas curtas com base no vento zonal máximo), herdado do ECMWF, que deveria mitigar parte dessas instabilidades. Recomenda-se verificar se esse mecanismo está ativo.

---

## 🧪 Testes Complementares

Um teste auxiliar com **inicialização `WarmStart` (initlz = -3)** também foi executado com sucesso em `DeltaT = 220 s` e `DeltaT = 200 s`, reforçando a eficácia do ajuste.

---

## 🚀 Procedimentos Recomendados para Pré-Operação

1. **Definir `DeltaT = 200 s`** na configuração do BAM/SMNA (T299L64).
2. **Executar a próxima simulação com `ColdStart`**, garantindo estado inicial limpo.
3. **Confirmar a frequência das chamadas de radiação:**

   * Se as chamadas ocorrerem em múltiplos de 1 h ou 3 h, assegurar que o `DeltaT` seja submúltiplo de 3600 s.
4. **Verificar ativação do controle de CFL** no código-fonte ou nas configurações dinâmicas.
5. **Monitorar a estabilidade** durante integrações longas (> 24 h) e ciclos sucessivos de assimilação.

---

## 🧾 Conclusão

O ajuste do **DeltaT para 200 s** apresentou o melhor resultado, eliminando o *segmentation fault* e estabilizando a integração, mesmo após várias horas de simulação.
Apesar do sucesso empírico, recomenda-se conduzir uma **análise mais profunda sobre a sensibilidade do timestep** e o comportamento do controle de CFL em diferentes condições atmosféricas.

---

**Autor:** João Gerd Zell de Mattos
**Colaborações:** Kubota, Bonatti
**Data:** Fevereiro de 2025

---
