# SMNA – Sistema de Modelagem Numérica e Assimilação

O **SMNA** é o acoplamento do Modelo Atmosférico Global Brasileiro (**BAM – Brazilian Global Atmospheric Model**) em sua versão com coordenada vertical híbrida com o Sistema de Assimilação de Dados **GSI – Global Statistical Interpolation**.

O SMNA compõe o sistema global de modelagem e assimilação de dados do **CPTEC/INPE**, integrando modelagem numérica atmosférica com assimilação variacional de observações.

A concepção inicial do sistema foi proposta no âmbito do desenvolvimento do acoplamento BAM + GSI. Ajustes e evoluções foram incorporados ao longo do desenvolvimento, conforme documentado nas versões do sistema.

---

## Componentes do SMNA

O sistema SMNA é composto pelos seguintes componentes principais:

- **BAM** – Brazilian Global Atmospheric Model  
- **GSI** – Gridpoint Statistical Interpolation  
- **SPCON** – Sistema de Previsão por Conjunto Global  

Cada componente possui arquitetura própria, mas o SMNA organiza o fluxo integrado de:

1. Preparação de dados
2. Assimilação variacional
3. Integração do modelo atmosférico
4. Pós-processamento
5. Ciclagem operacional

---

## Estrutura do Repositório

O repositório contém:

- Scripts de execução do ciclo de assimilação
- Configurações do GSI
- Fluxo de execução do BAM
- Integração com ambiente HPC (XC50/EGEON)
- Templates e arquivos de controle
- Controle de versões históricas do sistema

---

## Planejamento de Versões

Principais versões históricas do SMNA:

- **SMNA 2.2.0** – BAM híbrido + IBIS + assimilação de dados  
- **SMNA 2.3.0** – Inclusão de novo pré-processamento e correções no modelo  
- **SMNA 2.4.0** – Ajustes na minimização da função custo e correção de bias  
- **SMNA 2.5.0** – Introdução da matriz de covariâncias híbrida  
- **SMNA 2.6.0** – Consolidação da assimilação de superfície  
- **SMNA 2.7.0** – Ajustes adicionais na matriz de erro do BAM híbrido  
- **SMNA Oper** – Versão operacional do sistema  

Tags Git correspondentes estão disponíveis neste repositório.

---

## Ambiente de Execução

O SMNA é executado em ambiente HPC, com suporte a:

- Execução paralela MPI
- Ciclo variacional
- Execução operacional
- Ambiente XC50/EGEON

---

## Histórico

Este repositório é resultado da migração oficial do histórico SVN institucional do SMNA para Git, preservando:

- Histórico completo de commits
- Autores originais
- Branches históricos
- Tags de versões oficiais

---

## Organização

Grupo de Assimilação de Dados (GAD)  
Divisão de Modelagem Numérica do Sistema Terrestre (DIMNT)  
Centro de Previsão de Tempo e Estudos Climáticos (CPTEC/INPE)
