# Métodos de Controle de Qualidade (QC) no GSI

Baseado no manual do **GSI (Gridpoint Statistical Interpolation)**, os principais métodos de **controle de qualidade (QC)** aplicados aos dados observacionais no processo de assimilação incluem:

---

## 1. Caracterização e Atribuição de Erros de Observação
- Os erros associados às observações são continuamente ajustados e calibrados.  
- Inicialmente tratados como não correlacionados, mas avanços recentes no GSI incluem erros **espectralmente correlacionados**, usando estimativas baseadas em métodos de **Desroziers et al. (2005)**.

---

## 2. Controle de Qualidade Variacional (VarQC)
- Implementado para dados in situ.  
- Utiliza distribuições robustas como a **Huber norm**.  
- O VarQC não rejeita observações: **repondera** aquelas com grandes inovações (*obs – background*), reduzindo seu peso na análise.  
- Isso permite flexibilizar ou substituir o *gross check*, aproveitando mais dados sem comprometer a robustez.

---

## 3. Checagens Tradicionais de Consistência
Aplicadas em paralelo ao VarQC:

- **Gross Check** (*checagem grosseira*): remove observações com desvios muito grandes em relação ao background.  
- **Buddy Check**: compara observações vizinhas (espacial/temporal) para detecção de outliers.  
- **Thinning / Superobbing**: redução ou agregação de dados redundantes (especialmente satélite e radar).

---

## 4. Bias Correction (Correção de Viés)
- Fundamental para radiâncias de satélite.  
- Inclui **correção de viés variacional (VarBC)** combinada com QC.  
- Evita erros sistemáticos ao usar sensores orbitais.

---

## 5. Filtros Específicos por Tipo de Dado
Exemplos:

- **Radiossondas de alta resolução (BUFR):** ajustes específicos antes do uso operacional.  
- **GNSS-RO:** QC modernizado para lidar com fortes gradientes na PBL.  
- **Radar e AMVs:** filtros próprios, incluindo correções de erro correlacionado.

---

## Resumo
O **GSI** combina métodos clássicos (gross/buddy checks, filtros específicos), atribuição adaptativa de erros, correção de vieses e, mais recentemente, o **VarQC**, que permite aproveitar um conjunto maior de observações, mas com ponderação diferenciada segundo sua confiabilidade.

---

## Tabela Comparativa de Métodos de QC (Manual do GSI)

| **Método de QC**                         | **Descrição**                                                                                                                                                  | **Vantagens**                                                                                      | **Limitações**                                                                                  |
| ---------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| **VarQC (Quality Control Variacional)**  | Usa funções robustas (ex.: Huber norm) para reponderar observações com grandes inovações (obs – background). Não rejeita dados, apenas reduz o peso na análise. | Aproveita mais observações; robustez contra outliers; substitui em parte o *gross check*.           | Pode incluir dados de baixa qualidade, ainda que com peso reduzido.                             |
| **Gross Check**                          | Elimina observações com desvios muito grandes em relação ao campo de background.                                                                               | Simples e eficiente para remover outliers extremos.                                                | Pode descartar observações válidas em regiões de forte gradiente (frentes, convecção).          |
| **Buddy Check**                          | Compara observações vizinhas (espacial/temporal) para detectar inconsistências.                                                                                | Aumenta a confiabilidade ao rejeitar dados isolados incoerentes.                                   | Requer densidade suficiente de observações; pode falhar em regiões de baixa cobertura.          |
| **Thinning / Superobbing**               | Reduz ou agrega observações redundantes, especialmente de satélite e radar.                                                                                    | Diminui correlação entre erros; reduz custo computacional.                                         | Perda de parte da informação de alta resolução.                                                 |
| **Atribuição de Erros de Observação**    | Define e ajusta erros observacionais, incluindo métodos para erros correlacionados (Desroziers et al. 2005).                                                   | Melhor representação estatística; aumenta consistência da análise.                                 | Estimativas dependem de amostras representativas; erros podem ser mal caracterizados.           |
| **Bias Correction (Correção de Viés)**   | Ajusta sistematicamente observações (especialmente radiâncias) antes da assimilação. Pode ser variacional (VarBC).                                             | Essencial para uso de satélites; reduz erros sistemáticos.                                         | Requer monitoramento constante; risco de “aprender” erros do modelo.                            |
| **Filtros Específicos por Tipo de Dado** | Aplicações particulares: radiossondas (BUFR), GNSS-RO, AMVs, radar, etc., cada um com rotinas próprias de QC.                                                  | Ajusta o QC às características de cada sensor; aumenta a utilidade dos dados.                      | Complexidade maior; depende de conhecimento especializado de cada observação.                   |

---

# QC no GSI – detalhamento técnico com equações e referências

Abaixo está um “deep-dive” dos métodos de **controle de qualidade (QC)** mais usados no GSI, com ênfase na formulação variacional, onde a função-custo de DA é tipicamente

$$
J(\delta x)=\tfrac12\,\delta x^\mathrm{T}\mathbf{B}^{-1}\delta x\;+\;\tfrac12\,\sum_{i} \rho_i\!\left(r_i\right),
\quad
r_i=\frac{d_i-\mathbf{H}_i\,\delta x}{\sigma_i},
$$

onde $d_i=y_i-\mathbf{H}_i x_b$ é a inovação (obs–background), $\sigma_i$ é o desvio-padrão do erro de observação, $\mathbf{B}$ é a covariância de erro de background e $\rho_i(\cdot)$ define o “perfil” de QC (quadrático gaussiano, Huber, etc.). ([dtcenter.org][1])

---

## 0) Especificação do erro de observação $\mathbf{R}$ (não/ correlacionado)

### Caso padrão (diagonal)

No GSI, muitos tipos de dados usam $\mathbf{R}=\mathrm{diag}(\sigma_1^2,\dots,\sigma_n^2)$, com $\sigma_i$ lidos de tabelas (p.ex., `convinfo`, `satinfo`) e ajustes automáticos por qualidade, altitude, etc. Parâmetros de erro e *gross check* por tipo/subtipo são configuráveis em `convinfo` (p.ex., `cgross`, `cermax`, `cermin`). ([dtcenter.org][2])

### Erros correlacionados (radiâncias multicanais, *all-sky*)

Para satélites hiperespectrais (IASI, CrIS), há inter-correlações entre canais. Uma forma é escrever $\mathbf{R}=\mathbf{D}^{1/2}\mathbf{C}\mathbf{D}^{1/2}$, com $\mathbf{C}$ a matriz de correlação. A “branqueação” dos resíduos usa, por exemplo, fatoração de Cholesky $\mathbf{C}=\mathbf{L}\mathbf{L}^\mathrm{T}$ e resíduos transformados $\tilde r=\mathbf{L}^{-1}(\mathbf{D}^{-1/2}(d-\mathbf{H}\delta x))$. Estudos operacionais no Met Office/ECMWF mostraram ganhos ao diagnosticar $\mathbf{R}$ com correlações via estatísticas de $(O\!-\!B)$ e $(O\!-\!A)$. ([RMETS][3])

Uma ferramenta prática para aferir $\mathbf{R}$ é o diagnóstico de Desroziers et al. (2005), que fornece relações de consistência (quase “custo zero”) para $\mathbf{R}$, $\mathbf{B}$ e erro de análise em espaço de observação, a partir de $(O\!-\!B)$ e $(O\!-\!A)$. ([meteo.physik.lmu.de][4])

> Status no ecossistema NOAA/NCEP: diretrizes recentes enfatizam expandir o uso de VarQC (com **Huber**) e enfrentar erros correlacionados — sobretudo em condições *all-sky*. ([repository.library.noaa.gov][5])

---

## 1) VarQC (Quality Control Variacional)

### 1.1 Formulação “AJ99” (mistura Gaussiana + uniforme)

O VarQC clássico modela o erro de observação como mistura: com probabilidade $1-p$ o erro é Gaussiano ($\mathcal N(0,\sigma^2)$) e com probabilidade $p$ vem de uma distribuição “plana” (gross error). Isso induz uma penalização não-quadrática equivalente a **reponderar** cada observação por um peso $w_i\in(0,1]$:

$$
w_i \;=\;
\frac{(1-p)\,\phi(r_i)}{(1-p)\,\phi(r_i)+p\,c},\qquad
J_o=\tfrac12\sum_i w_i\,r_i^2,
$$

onde $\phi$ é a densidade Normal padrão e $c$ é uma constante da “cauda plana”. Pesos pequenos (inovações grandes) reduzem a influência da observação sem descartá-la. ([RMETS][6])

### 1.2 Norma de **Huber** (implementação moderna)

Na prática operacional europeia e mais recentemente no ecossistema NOAA, usa-se a **Huber norm**:

$$
\rho_\kappa(r)=
\begin{cases}
\tfrac12\,r^2,& |r|\le \kappa,\\[2pt]
\kappa\,|r|-\tfrac12\,\kappa^2,& |r|>\kappa,
\end{cases}
\qquad
w(r)=\frac{\partial^2 \rho}{\partial r^2}=
\begin{cases}
1,& |r|\le \kappa,\\
0,& |r|>\kappa \text{ (no sentido de “linear”)},
\end{cases}
$$

que é **quadrática** na região central (gaussiana) e **linear** nas caudas (robusta). Impactos positivos para dados *in situ* foram documentados no ECMWF; a NOAA/NCEP tem avançado com VQC baseado em Huber (NCEP-VQC). ([ECMWF][7])

> No GSI, VarQC com Huber foi incorporado para *in situ* (e segue em expansão), mantendo o espírito de **reponderação** (em vez de rejeição “dura”). ([repository.library.noaa.gov][5])

---

## 2) *Gross Check* (checagem grosseira)

Objetivo: rejeitar valores com $|d_i|$ muito grandes face à incerteza combinada. Um critério típico é

$$
\left|\frac{d_i}{\sqrt{S_i + \sigma_i^2}}\right| \; \le\; \texttt{gross},
$$

em que $S_i$ é a variância do *background* (ou *ensemble spread* em filtros de conjunto). No GSI, os limiares por tipo/subtipo residem em `convinfo` (campos `cgross`, `cermax`, `cermin`). Para *Ensemble* (EnSRF/EnKF), há parâmetros equivalentes (p.ex., `sprd_tol`, `varqc`, `huber`). ([dtcenter.org][2])

---

## 3) *Buddy Check* (consistência espaço-temporal)

Verifica a consistência local comparando uma observação $i$ com “vizinhos” $j$ em janelas espaço-temporais:

$$
T_{ij}=\frac{(d_i-d_j)}{\sqrt{\sigma_i^2+\sigma_j^2 + 2\,\mathrm{Cov}(e_i,e_j)}},
\quad \text{rejeita se}\; |T_{ij}|>\tau.
$$

A formulação moderna é frequentemente Bayesiana/multivariada (Ingleby & Lorenc), combinando informação de vários vizinhos para inferir a probabilidade de “erro grosso”. ([RMETS][8])

---

## 4) *Thinning* e **Superobbing**

**Thinning** reduz correlação de erro e custo computacional impondo malha mínima (e.g., `ithin_conv`, `rmesh_conv` em `convinfo`). **Superobbing** agrega $N$ observações redundantes numa célula e atualiza o erro agregado:

$$
\mathbf{Var}(\bar e)=\frac{1}{N^2}\sum_{i,j}\mathbf{Cov}(e_i,e_j)
\;\approx\;
\frac{\sigma^2}{N}\,\bigl[\,1+(N-1)\rho\,\bigr],
$$

se $\sigma^2$ e a correlação média $\rho$ forem homogêneas. Assim, se $\rho>0$, o ganho é **menor** que $1/\sqrt{N}$, justificando *thinning/superobbing* com atenção à correlação. ([Tellus A][9])

---

## 5) Correção de Viés (**VarBC**) – radiâncias

A correção de viés é crítica para radiâncias. No **VarBC**, o viés é modelado como $b=\mathbf{Z}\beta$ (preditores: constante, ângulo de varredura, termos de massa de ar, etc.), e $\beta$ entra como **variável de controle** adicional:

$$
J(\delta x,\beta)=\tfrac12\,\delta x^\mathrm{T}\mathbf{B}^{-1}\delta x
+\tfrac12\sum_i \frac{\bigl(y_i-\mathbf{H}_i(x_b+\delta x)-\mathbf{z}_i^\mathrm{T}\beta\bigr)^2}{\sigma_i^2}
+\tfrac12\,(\beta-\beta_b)^\mathrm{T}\mathbf{B}_\beta^{-1}(\beta-\beta_b).
$$

Isso permite ajuste adaptativo e “ancoragem” por observações estáveis (superfície, radiossondas) para evitar que a correção aprenda viés do modelo. ([leg.ufpr.br][10])

---

## 6) Filtros específicos por tipo de dado (exemplos)

### 6.1 Radiossondas (*BUFR* de alta resolução)

Checagens de consistência vertical (p.ex., razão de variação térmica), limites físicos (UR, vento), janelas temporais e de *tracking* da estação; no GSI, a aceitação e as janelas de tempo/altura são tratadas nos leitores/rotinas de pré-QC e no próprio *gross check*. ([dtcenter.org][1])

### 6.2 GNSS-RO (bending angle/refratividade)

QC fortemente dependente da altitude (camada limite) e de indicadores de qualidade de *retrieval* (p.ex., **LSW** – *logarithm of spectral width*). Esquemas **LSW-dependentes** rejeitam perfis com viés sistêmico na baixa troposfera. A assimilação direta de *bending angle* com checagens à la $\chi^2$ é padrão. ([RMETS][11])

### 6.3 AMVs (vetores de movimento atmosférico)

Uso de **Quality Indicator (QI)**/Expected Error (EE) e filtros de altura/consistência dinâmica (coerência do campo, *speed bias*). Limiarizações de QI (p.ex., QI ≥ 60) e *thinning* são práticas correntes. ([journals.ametsoc.org][12])

### 6.4 Radar (ventos radiais/reflectividade)

Checagens de *dealiasing*, eco de solo/ruído, filtros por intensidade/altura e **superobbing** polar para reduzir correlações. A especificação de $\mathbf{R}$ e *gross check* depende do tipo (radial wind vs Z), com forte ênfase em representatividade (hidrometeoros/convectivo). (Ver diretrizes *all-sky* para correlações situação-dependentes.) ([repository.library.noaa.gov][5])

---

## 7) Itens práticos no **GSI**

* **Parâmetros em `convinfo`:** `cgross` (limiar *gross*), `cermax/cermin` (faixas de erro), `cvar_b/cvar_pg` (VarQC), `ithin_conv`, `rmesh_conv` (thinning). ([nco.ncep.noaa.gov][13])
* **Listas de uso/rejeição:** `uselist`/`rejection list` controlam habilitação por estação/sensor além dos QC internos. ([dtcenter.org][14])
* **Diagnósticos para ajuste de $\mathbf{R}$:** use estatísticas $(O\!-\!B)$ e $(O\!-\!A)$ (método de **Desroziers**) para checar consistência entre pesos/erros. ([meteo.physik.lmu.de][4])

---
Aqui está um **cheat-sheet em Markdown** consolidando os métodos de QC do **GSI** (baseado no manual e no documento NOAA/NCEP/EMC Strategy), com exemplos de parâmetros típicos em `convinfo` e blocos de `namelist` relevantes para VarQC/Huber, listas de rejeição e *thinning*.

---

# 📑 GSI QC Cheat-Sheet

## 🔹 Métodos de QC no GSI

* **VarQC/Huber norm**: repondera inovações extremas sem rejeição dura.
* **Gross check**: rejeição com limiar de desvio normalizado.
* **Buddy check**: consistência espaço-temporal entre vizinhos.
* **Thinning / Superobbing**: redução de dados redundantes.
* **Bias correction (VarBC)**: correção de viés em radiâncias.
* **Filtros específicos**: AMVs (QI), GNSS-RO (LSW), radiossondas BUFR, radar.

---

## 🔹 Exemplo de Linha `convinfo`

```text
! varname  type  kx   error   gross   cvar_b   cvar_pg   ithin_conv   rmesh_conv
  ps       120   120  1.50    10.0    1.0      1.0       1            145.0
  t        120   120  1.00    7.0     1.0      1.0       1            145.0
  uv       220   220  1.50    8.0     1.0      1.0       1            145.0
```

* **error**: σ da observação.
* **gross**: limiar do *gross check* (|O-B|/σ).
* **cvar\_b / cvar\_pg**: parâmetros VarQC.
* **ithin\_conv / rmesh\_conv**: opções de *thinning* (ativação e distância em km).

---

## 🔹 `uselist` e `rejection list`

Arquivos externos listam estações/sensores a usar ou rejeitar:

```text
! uselist.txt
  72451
  72452

! rejection_list.txt
  72295
  72297
```

* **uselist** → força inclusão (se passar QC).
* **rejection** → rejeita independentemente do QC.

---

## 🔹 `namelist` – VarQC (Huber)

```fortran
&qc_var_settings
  varqc      = .true.      ! ativa VarQC
  huber      = .true.      ! usa norma de Huber
  huber_k    = 2.0         ! parâmetro de corte (|O-B|/σ)
  gross_fac  = 10.0        ! fator do gross check
  sprd_tol   = 5.0         ! tolerância para ensemble spread
/
```

* **huber\_k** controla a transição entre quadrático e linear.
* Valores típicos: 1.5 – 2.5 (ajustar por tipo de dado).

---

## 🔹 `namelist` – Thinning/Superobbing

```fortran
&thinning_settings
  ithin_conv = 1        ! ativa thinning conv.
  rmesh_conv = 145.0    ! malha mínima (km)
  ithin_rad  = 1        ! thinning em radiâncias
  rmesh_rad  = 75.0
  superob    = .true.   ! ativa superobbing
/
```

* Distâncias típicas: **145 km** (convencional), **75 km** (radiância).
* *Superobbing* recomendado para radar e satélite de alta resolução.

---

## 🔹 Fluxo Prático

1. Ajuste **σ e gross** em `convinfo`.
2. Ative **VarQC/Huber** no `namelist`.
3. Defina listas `uselist`/`rejection` conforme monitoramento.
4. Ajuste *thinning* por tipo de dado (convencional, satélite, radar).
5. Monitore estatísticas (O-B, O-A) e use diagnósticos de Desroziers para validar erros.

---

## Referências selecionadas (porta de entrada)

* **VarQC e Huber:** Andersson & Järvinen (1999); Tavolato & Isaksen (2014/2015). ([RMETS][6])
* **Diagnóstico de erros (R, B):** Desroziers et al. (2005). ([Sistema de Dados de Astrofísica][15])
* **Erros correlacionados em radiâncias:** Stewart et al. (2014); Bormann et al. (2015). ([RMETS][16])
* **GSI – Guias do Usuário (v3.5/v3.6) e *Advanced Guide*:** DTCenter (2016–2017). ([dtcenter.org][1])
* **Diretrizes NOAA/NCEP recentes:** *Data Assimilation Strategy for NOAA/NWS/NCEP/EMC* (2024). ([repository.library.noaa.gov][5])
* **Buddy/Bayesian QC:** Ingleby & Lorenc (1993). ([RMETS][8])
* **GNSS-RO QC (LSW, bending angle):** Healy & Thépaut (2006); Liu et al. (2018). ([RMETS][11])
* **AMVs (QI/EE):** Bedka et al. (2005); Wanzong et al. (2012). ([journals.ametsoc.org][12])

---

## Dicas finais de uso

1. **Comece simples**: *gross check* e erros $\sigma$ coerentes com as estatísticas $(O\!-\!B)$. Valide com Desroziers. ([meteo.physik.lmu.de][4])
2. **Habilite VarQC (Huber)** para *in situ* quando houver caudas pesadas; ajuste $\kappa$ (*tuning* conservador). ([ECMWF][7])
3. **Radiâncias**: cuide primeiro de **VarBC**; só então explore $\mathbf{R}$ correlacionado (comece por submatrizes/channels-clusters para evitar condicionamento ruim). ([leg.ufpr.br][10])
4. **Dados densos** (AMV, radar): *thinning/superobbing* com célula e erro ajustados à correlação esperada ($\rho$). ([Tellus A][9])


[1]: https://dtcenter.org/sites/default/files/community-code/gsi/docs/users-guide/GSIUserGuide_v3.5.pdf?utm_source=chatgpt.com "Users Guide"
[2]: https://dtcenter.org/sites/default/files/GSIUserGuide_v3.6.pdf?utm_source=chatgpt.com "User's Guide Version 3.6"
[3]: https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.2211?utm_source=chatgpt.com "Estimating interchannel observation‐error correlations for IASI ..."
[4]: https://www.meteo.physik.lmu.de/DokuWiki/lib/exe/fetch.php?media=lscraig%3Aherz%3Al95%3Alectures%3Adesroziers-qjrms131-2005_en.pdf&utm_source=chatgpt.com "Diagnosis of observation, background and analysis-error ..."
[5]: https://repository.library.noaa.gov/view/noaa/56694/noaa_56694_DS1.pdf?utm_source=chatgpt.com "Data Assimilation Strategy for NOAA/NWS/NCEP/EMC"
[6]: https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1002/qj.49712555416?utm_source=chatgpt.com "Variational quality control"
[7]: https://www.ecmwf.int/sites/default/files/elibrary/2014/12572-use-huber-norm-observation-quality-control-ecmwf-4d-var.pdf?utm_source=chatgpt.com "On the use of a Huber norm for observation quality control ..."
[8]: https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.49711951316?utm_source=chatgpt.com "Bayesian quality control using multivariate normal distributions"
[9]: https://a.tellusjournals.se/articles/10.3402/tellusa.v66.23294?utm_source=chatgpt.com "Estimating correlated observation error statistics using an ..."
[10]: https://www.leg.ufpr.br/~eder/Artigos/Bias/Auligne_2007.pdf?utm_source=chatgpt.com "Adaptive bias correction for satellite data in a numerical ..."
[11]: https://rmets.onlinelibrary.wiley.com/doi/10.1256/qj.04.182?utm_source=chatgpt.com "Assimilation experiments with CHAMP GPS radio occultation ..."
[12]: https://journals.ametsoc.org/view/journals/apme/44/11/jam2264.1.xml?utm_source=chatgpt.com "Application of Satellite-Derived Atmospheric Motion Vectors ..."
[13]: https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/rap.v5.1.20/sorc/rap_gsi.fd/src/gsi/convinfo.f90?utm_source=chatgpt.com "convinfo.f90"
[14]: https://dtcenter.org/sites/default/files/community-code/gsi/docs/users-guide/AdvancedGSIUserGuide_v3.5.0.0.pdf?utm_source=chatgpt.com "Advanced User's Guide"
[15]: https://ui.adsabs.harvard.edu/abs/2005QJRMS.131.3385D/abstract?utm_source=chatgpt.com "Diagnosis of observation, background and analysis-error ..."
[16]: https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1002/qj.2211?utm_source=chatgpt.com "Estimating interchannel observation-error correlations for IASI ..."
