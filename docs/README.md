# Documentação Técnica — SMNA

Este diretório reúne documentação técnica consolidada do SMNA (Sistema de Modelagem Numérica e Assimilação), incluindo aspectos de:

* Configuração do BAM
* Parametrizações do GSI
* Controle de qualidade
* Execução operacional
* Compilação
* Ambientes HPC

---

# 1. Modelo BAM

## Ajuste do DeltaT

`BAM_DELTAT.md`
Análise detalhada da instabilidade numérica no BAM (T299L64), diagnóstico do erro de segmentação e testes com diferentes valores de `DeltaT`, incluindo implicações na condição CFL e chamadas de radiação.

---

# 2. Sistema de Assimilação GSI

## Parâmetros do OBS_INPUT

`OBS_INPUT.md`
Análise técnica dos parâmetros `dval`, `dthin` e `dsfcalc`, baseada em inspeção direta do código-fonte do GSI.

---

## Controle de Qualidade no GSI

`QC_GSI.md`
Descrição dos métodos de controle de qualidade no GSI: VarQC, Huber norm, gross check, buddy check, thinning e bias correction.

`qc_marker.md`
Documentação das marcas de controle de qualidade (`idqc`, `iuse`, `lim_qm`, `noiqc`) no fluxo PREPBUFR → GSI.

---

## PCGSOI – Cálculo do Stepsize

`pcgsoi.md`
Descrição matemática e estrutural do algoritmo `stpcalc` utilizado no método PCGSOI no GSI, incluindo representação interna dos coeficientes e estrutura diferencial da matriz `pbc`.

---

# 3. PREPBUFR

## Leitura de Mensagens

`Read_Prepbufr_Messages.md`
Lista e interpretação das mensagens do código `read_prepbufr.f90`, incluindo erros fatais, avisos e mensagens informativas.

---

# 4. Implementação Pré-Operacional

## SMNA no Egeon

`SMNA_egeon.md`
Guia para obtenção, configuração, compilação e execução do SMNA no cluster Egeon.

## Implementação Pré-Operacional Geral

`SMNA_pre-operacional.md`
Estratégia de implantação pré-operacional no XC50 com preparação paralela para o Egeon.

## Tag Estável SMNA_v3.0.1

`SMNA_v3.0.1.md`
Guia específico para a tag `SMNA_v3.0.1`, com início rápido, bootstrap do ciclo e execução com spin-up de bias.

---

# 5. Compilação e Ferramentas

## PostAlt

`PostAlt_comp.md`
Passo-a-passo de compilação do PostAlt (local, Egeon e XC50), incluindo dependências `sharedLibs`, Autotools e diferenças de ambiente.

## Diretórios entre Clusters

`paths.md`
Correspondência de diretórios entre BASTOS, EGEON e XC50, com recomendações para scripts portáveis.

---

# 6. Bias Correction

## Guia de Spin-up de Satbias

`Bc.md`
Guia para spin-up de coeficientes de correção de viés (`satbias` e `satbias_angle`) no GSI, incluindo checagens de consistência e script de automação.
