# Implementação Pré-Operacional do SMNA — Egeon

> **Resumo:** Este guia explica, de ponta a ponta, como **obter**, **configurar**, **compilar** e **executar** o ciclo pré-operacional do **SMNA** (BAM + GSI) no **Egeon**.

## Sumário

* [1. Visão geral](#1-visão-geral)
* [2. Escopo e estratégia](#2-escopo-e-estratégia)
* [3. Obtenção do sistema (SVN)](#3-obtenção-do-sistema-svn)

  * [Passo a passo (Egeon)](#passo-a-passo-egeon)
  * [Estrutura esperada após o download](#estrutura-esperada-após-o-download)
  * [Comparação entre métodos](#comparação-entre-métodos)
* [4. Instalação e configuração](#4-instalação-e-configuração)

  * [4.1 Editar o arquivo de paths](#41-editar-o-arquivo-de-paths)
  * [4.2 Configurar ambiente inicial](#42-configurar-ambiente-inicial)
  * [4.3 Compilar os componentes](#43-compilar-os-componentes)
  * [4.4 Copiar arquivos fixos](#44-copiar-arquivos-fixos)
  * [4.5 Verificar executáveis](#45-verificar-executáveis)
  * [4.6 Rodar caso de teste](#46-rodar-caso-de-teste)
* [5. Compilação (Intel no Egeon)](#5-compilação-intel-no-egeon)

  * [Carregar módulos](#carregar-módulos)
  * [Compilar](#compilar)
  * [Verificar executáveis gerados](#verificar-executáveis-gerados)
* [6. Execução: testcase e ciclo](#6-execução-testcase-e-ciclo)

  * [6.1 Preparação](#61-preparação)
  * [6.2 Testcase](#62-testcase)
  * [6.3 Pré-processamento do BAM (uma única vez por resolução)](#63-pré-processamento-do-bam-uma-única-vez-por-resolução)

    * [`runPre` — referência rápida de opções](#631-runpre--referência-rápida-de-opções)
  * [6.4 First Guess (9 h) e restarts](#64-first-guess-9-h-e-restarts)

    * [`runModel` — referência rápida](#641-runmodel--referência-rápida)
  * [6.5 Ciclo de assimilação (operacional)](#65-ciclo-de-assimilação-operacional)
  * [6.6 Ciclo real (outra data): entradas NCEP/GDAS/SST](#66-ciclo-real-outra-data-entradas-ncepgdassst)

    * [Exemplo simples](#exemplo-simples)
    * [Opções disponíveis](#opções-disponíveis)
    * [Exemplos práticos](#exemplos-práticos)
* [7. Resultados e saídas geradas](#7-resultados-e-saídas-geradas)

  * [Análises](#análises)
  * [Restarts](#restarts)
  * [First Guess](#first-guess)
  * [Previsões (GRIB)](#previsões-grib)
  * [Tabela-resumo de saídas (Egeon)](#tabela-resumo-de-saídas-no-egeon-genérico)
  * [Exemplo concreto — Testcase 2019](#exemplo-concreto--testcase-2019-ciclo2019111500)
* [8. Dicas e solução de problemas](#8-dicas-e-solução-de-problemas)
* [Anexo: Início rápido (cheatsheet)](#anexo-início-rápido-cheatsheet)

---

## 1. Visão geral

O **SMNA** integra:

* **BAM** (*Brazilian Atmospheric Model*): previsão numérica
* **GSI** (*Gridpoint Statistical Interpolation*): assimilação de dados

Para o ciclo, são necessários:

1. **Fluxo de dados**, 2) **First guess** (estimativa *a priori*), 3) **Estatísticas de erro**.
   Plano atual: dados NCEP (fase 1), BAM com janela de **9 h**, e GSI pré-operacional.

---

## 2. Escopo e estratégia

Alvo: **Egeon**. Scripts (`config_smg.ksh`, `smg_setup.sh`) ajustados para este ambiente.

---

## 3. Obtenção do sistema (SVN)

Código-fonte do **SMNA/SMG** no **SVN do CPTEC**:

* **`tag`** → versões **estáveis** (usuários)
* **`trunk`** → desenvolvimento principal (contribuidores)
* **`branch`** → linhas de trabalho específicas

### Passo a passo (Egeon)

1. **SSH no Egeon**

```bash
ssh usuario@egeon.cptec.inpe.br
```

2. **Diretório HOME**

```bash
cd /home/${USER}
```

3. **Export de uma tag (ex.: `SMNA_v2.3.1`)**

```bash
svn export https://svn.cptec.inpe.br/smna/tag/SMNA_v2.3.1
```

*(Sem vínculo de versionamento; indicado para uso final.)*

4. **Listar tags**

```bash
svn list https://svn.cptec.inpe.br/smna/tag/
```

5. **Checkout do trunk (desenvolvimento)**

```bash
svn co https://svn.cptec.inpe.br/smna/trunk/SMNA
```

6. **Checkout de um branch (ex.: `SMNA_v3.0.0.t11889`)**

```bash
svn co https://svn.cptec.inpe.br/smna/branch/SMNA_v3.0.0.t11889
```

### Estrutura esperada após o download

* **Tag (export):** `/home/${USER}/SMNA_v2.3.1/SMG`
* **Trunk (checkout):** `/home/${USER}/SMNA/SMG`
* **Branch (checkout):** `/home/${USER}/SMNA_v3.0.0.t11889/SMG`

### Comparação entre métodos

| Método                | Comando exemplo                                                   | Indicação de uso                                               |
| --------------------- | ----------------------------------------------------------------- | -------------------------------------------------------------- |
| **Export (tag)**      | `svn export https://svn.cptec.inpe.br/smna/tag/SMNA_v2.3.1`       | Usuários finais: usar versão estável, sem versionamento local. |
| **Checkout (trunk)**  | `svn co https://svn.cptec.inpe.br/smna/trunk/SMNA`                | Contribuidores no desenvolvimento principal.                   |
| **Checkout (branch)** | `svn co https://svn.cptec.inpe.br/smna/branch/SMNA_v3.0.0.t11889` | Trabalhar/acompanhar uma linha específica.                     |

> **Nota importante**
>
> * Usuários: prefira **tags** (export).
> * Devs: use **checkout** (trunk/branch).
> * Este guia usa o **branch** `SMNA_v3.0.0.t11889` nos exemplos.

---

## 4. Instalação e configuração

### 4.1 Editar o arquivo de paths

Entre na pasta do **SMG** e edite `egeon_paths.conf`:

```bash
cd /home/${USER}/<NOME_SMNA>/SMG
nano etc/mach/egeon_paths.conf
```

Ajuste `nome_smg` conforme o diretório baixado:

```
# Exemplo 1:
nome_smg    SMNA_v3.0.0.t11889/SMG

# Exemplo 2:
nome_smg    SMNA_TESTE/SMG
```

### 4.2 Configurar ambiente inicial

Cria a estrutura de diretórios e faz cópias/links simbólicos:

```bash
./config_smg.ksh configure
```

> **Nota**
>
> * Todas as funções aceitam `-v` (verbose) e `-q` (quiet).
> * Use `-y` para responder “yes” automaticamente às perguntas (execução não interativa).

### 4.3 Compilar os componentes

Comando básico:

```bash
./config_smg.ksh compile
```

**Verbosidade**: `-v`, `-q`
**Seleção**: `--all | --none | --gsi | --no-gsi | --bam | --no-bam | --inctime | --no-inctime | --ang | --no-ang`
**Outros**: `-y` (auto-yes), `--dry-run`, `--restore`, `-f/--fix`, `-c/--cmake X` ou `--cmake=X`

**Exemplos**

```bash
# Compilar tudo em modo verboso:
./config_smg.ksh compile -v --all

# Apenas GSI e BAM:
./config_smg.ksh compile --gsi --bam --no-ang --no-inctime

# Dry-run (sem compilar de fato):
./config_smg.ksh compile --all --dry-run
```

> **Nota (ANGUPDATE)**
> A correção de viés é feita **dentro do GSI** nesta versão; **ANGUPDATE não é necessário**.
> Recomenda-se:

```bash
./config_smg.ksh compile --all --no-ang
```

### 4.4 Copiar arquivos fixos

```bash
./config_smg.ksh copy_fixed_files
```

### 4.5 Verificar executáveis

```bash
./config_smg.ksh verify_executables
```

Aceita as mesmas flags de verbosidade/seleção usadas em `compile`.

### 4.6 Rodar caso de teste

```bash
./config_smg.ksh testcase -v
```

---

## 5. Compilação (Intel no Egeon)

### Carregar módulos

```bash
module purge
module load intel/2022.0.2
module load mpi/2021.5.1

export NETCDF_DIR='/mnt/beegfs/lib_intel/netcdf'
export NETCDF='/mnt/beegfs/lib_intel/netcdf'
export HDF_DIR='/mnt/beegfs/lib_intel/hdf5-1.12.1'
export PNETCDF='/mnt/beegfs/lib_intel/pnetcdf-1.12.3/PnetCDF'
export PIO='/mnt/beegfs/lib_intel/pio-2.5.4'
export PATH=${PATH}:/mnt/beegfs/lib_intel/netcdf/bin
export AEC_LIBRARY='/mnt/beegfs/lib_intel/libaec-v0.3.2/src'
export AEC_INCLUDE_DIR='/mnt/beegfs/lib_intel/libaec-v0.3.2/include'
```

### Compilar

```bash
cd /home/${USER}/SMNA_v3.0.0.t11889/SMG

# Se ainda não fez:
./config_smg.ksh configure -v -y

# Compilar e salvar log:
./config_smg.ksh compile -v 2>&1 | tee compile_smg.log

# Acompanhar:
tail -f compile_smg.log
```

### Verificar executáveis gerados

Executáveis são disponibilizados no **HOME** (pastas `build/`) e no **SUBMIT\_HOME** (execuções reais). Caminhos típicos:

```bash
# GSI
/home/${USER}/SMNA_v3.0.0.t11889/SMG/cptec/gsi/build/src/gsi/gsi.x
/mnt/beegfs/${USER}/SMNA_v3.0.0.t11889/SMG/cptec/bin/gsi.x
# (ANGUPDATE desnecessário nesta versão)

# IncTime
/mnt/beegfs/${USER}/SMNA_v3.0.0.t11889/SMG/cptec/bin/inctime

# Pré-processamento (BAM)
/home/${USER}/SMNA_v3.0.0.t11889/SMG/cptec/bam/pre/build/ParPre_MPI
/mnt/beegfs/${USER}/SMNA_v3.0.0.t11889/SMG/bam/pre/exec/ParPre_MPI

# BAM (modelo)
/home/${USER}/SMNA_v3.0.0.t11889/SMG/cptec/bam/model/build/ParModel_MPI
/mnt/beegfs/${USER}/SMNA_v3.0.0.t11889/SMG/bam/model/exec/ParModel_MPI

# Pós-processamento (BAM)
/mnt/beegfs/${USER}/SMNA_v3.0.0.t11889/SMG/bam/pos/exec/bam/pos/PostGrib
```

> **Notas**
>
> * `build/` (HOME) é útil para depuração.
> * Use os executáveis do **SUBMIT\_HOME** nas execuções reais.
> * Se algo faltou, verifique `compile_smg.log` e repita.

---

## 6. Execução: testcase e ciclo

Após configurar/compilar, prepare as entradas:

* **Testcase** (dados prontos, ex.: 2019)
* **Ciclo real** (copia GDAS/SST de uma data específica)

### 6.1 Preparação

```bash
./config_smg.ksh copy_fixed_files
```

> Durante o `configure` ocorre a **primeira cópia** dos fixos.
> `copy_fixed_files` usa **rsync** (rápido, só copia o que mudou).

### 6.2 Testcase

```bash
./config_smg.ksh testcase -v
```

Solicita o ano (ex.: `2019`) e copia para `${SUBMIT_HOME}/datainout/`.
*Não executa o ciclo* — use `runPre`, `runModel` ou `run_cycle.sh`.

> **Nota:** o *testcase* de **2019** usa o ciclo **2019111500** (exemplos abaixo seguem essa data).

### 6.3 Pré-processamento do BAM (uma única vez por resolução)

```bash
cd /home/${USER}/SMNA_v3.0.0.t11889/SMG/cptec/bam/run
./runPre -v -t 299 -l 64 -I 2019111500 -n 0 -p SMT -s -O -T -G -Gt Grid
```

**Função:** prepara campos fixos/grade (orografia, máscaras etc.).
Em geral, **uma vez por resolução**.

> **Importante:** dados atuais do **NCEP** já vêm em NetCDF com prefixo **gdas**.
> Esses são os **padrões** do `runPre`; não precisa informar `-Gp gdas -Gt Grid` salvo exceções.

#### 6.3.1 `runPre` — referência rápida de opções

* `-t <TRUNC>`, `-l <LEVELS>` — resolução
* `-I <YYYYMMDDHH>` — data base para organização
* `-n <PART>` — partição (use 0 se não aplicável)
* `-p <SUITE>` — física (ex.: `SMT`, `CPT`)
* `-s` `-O` `-T` `-G` — ativa etapas (solo, orog., tsfc, grades)
* `-Gp <GRIDP>` (padrão: `gdas`)
* `-Gt <GRIDT>` (padrão: `netcdf`)
* `-v` / `-h`

### 6.4 First Guess (9 h) e restarts

Gera **FG** e **restarts** do BAM:

```bash
cd /home/${USER}/SMNA_v3.0.0.t11889/SMG/cptec/bam/run
./runModel -das -np 128 -t 299 -l 64 -I 2019111500 -W 2019111509 -F 2019111509 -ts 3 -r -tr 6 -i 2 -p SMT -s sstwkl
```

#### 6.4.1 `runModel` — referência rápida

* `-das` — modo acoplado à assimilação
* `-np <MPI>` — n° processos MPI
* `-t/-l` — resolução
* `-I` / `-W` / `-F` — datas (início / fim warm-up / fim)
* `-ts <H>` — frequência das saídas
* `-r` — grava restarts; `-tr <H>` — espaçamento dos restarts
* `-i <ITER>` — iterações internas
* `-p <SUITE>` — física; `-s <TAG>` — forçantes/superfície
* `-v` / `-h`

### 6.5 Ciclo de assimilação (operacional)

```bash
cd /home/${USER}/SMNA_v3.0.0.t11889/SMG/run

# ajuda
./run_cycle.sh -h

# execução típica
./run_cycle.sh -t 299 -l 64 -gt 254 -p CPT -I 2019111500 -F 2019112018

# execução reconectável
nohup ./run_cycle.sh -t 299 -l 64 -gt 254 -p CPT -I 2019111500 -F 2019112018 \
  | sed $'s/\e\[[0-9;:]*[a-zA-Z]//g' > run_cycle.out &
tail -f run_cycle.out
```

> **Notas**
> • Pós-processamento GRIB do BAM vem desativado; habilite no `run_cycle.sh` (troque `"No"` por `"Yes"` na chamada final ao `run_model.sh`).
> • Para produtos **6 h** (em vez de 3 h), altere `kpds=10 → kpds=11` no `POSTIN-GRIB` (pasta `run`).
> • Garanta entradas GDAS/SST (use `copy_ncep_input` antes).

### 6.6 Ciclo real (outra data): entradas NCEP/GDAS/SST

#### Exemplo simples

```bash
./config_smg.ksh copy_ncep_input -v 2025091400
```

#### Opções disponíveis

* `--dry-run` — só valida/lista, sem copiar
* `-v/--verbose` — logs detalhados
* `--src-root DIR` — sobrescreve origem padrão (`/oper/dados/ioper/tempo/NCEP/input`)
* `--layout LAYOUT` — `ncep-gfs` (padrão) ou `pre-year` (árvore antiga por ano)

#### Exemplos práticos

```bash
# Ciclo padrão (layout ncep-gfs)
./config_smg.ksh copy_ncep_input 2013010100

# Validar sem copiar
./config_smg.ksh copy_ncep_input --dry-run 2024021000

# Usar fonte alternativa
./config_smg.ksh copy_ncep_input --src-root /mnt/beegfs/backup/NCEP/input 2024010100

# Layout antigo (pré-2019)
./config_smg.ksh copy_ncep_input --layout pre-year \
  --src-root "${public_bam}/PRE/datain" 2019010100
```

**Destino das cópias:**

```
${SMG}/datainout/bam/pre/datain
```

---

## 7. Resultados e saídas geradas

### Análises

* **Dir:** `/mnt/beegfs/${USER}/SMNA_v3.0.0.t11889/SMG/datainout/gsi/dataout/<DATA>`
* **Arquivo:** `GANLCPT<AAAAMMDD><HH>S.unf.TQ0299L064`

### Restarts

* **Dir:** `/mnt/beegfs/${USER}/SMNA_v3.0.0.t11889/SMG/datainout/bam/pos/dataout/<RESOLUCAO>/<DATA>`
* **Arquivos típicos:** `GFCTCPT...convclP<rank>`, `...outattP<rank>`, `...outmdtP<rank>`, `...sibprgP<rank>`

### First Guess

* **Dir:** idem Restarts
* **Arquivos:** `GANLCPT<ana><prev>F.dir.<res_mcga>`, `GANLCPT<ana><prev>F.fct.<res_mcga>`

### Previsões (GRIB)

* **Dir:** idem Restarts
* **Arquivos:** `GPOSCPT<ana><ana>P.icn.<res_mcga>.grb`, `...inz...`, `GPOSCPT<ana><prev>P.fct.<res_mcga>.grb`

> **Prévia em grade (full):** `GANLNMC<ana>S.unf.<res_mcga>.GrADS` em
> `${SUBMIT_HOME}/SMG/datainout/bam/pre/dataout`.

### Tabela-resumo de saídas no Egeon (genérico)

| Tipo                 | Diretório base                                              | Padrão de arquivos                                                 |
| -------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------ |
| **Análises**         | `${SUBMIT_HOME}/SMG/datainout/gsi/dataout/<DATA>`           | `GANLCPT<AAAAMMDD><HH>S.unf.TQ0299L064`                            |
| **Restarts**         | `${SUBMIT_HOME}/SMG/datainout/bam/pos/dataout/<RES>/<DATA>` | `GFCTCPT...F.unf.<res>.(convclP\|outattP\|outmdtP\|sibprgP)<rank>` |
| **First Guess**      | `${SUBMIT_HOME}/SMG/datainout/bam/pos/dataout/<RES>/<DATA>` | `GANLCPT<ana><prev>F.dir.<res_mcga>`, `...fct.<res_mcga>`          |
| **Previsões (GRIB)** | `${SUBMIT_HOME}/SMG/datainout/bam/pos/dataout/<RES>/<DATA>` | `GPOSCPT...icn/inz/fct.<res_mcga>.grb`                             |
| **Prévia (GrADS)**   | `${SUBMIT_HOME}/SMG/datainout/bam/pre/dataout/`             | `GANLNMC<ana>S.unf.<res_mcga>.GrADS`                               |

### Exemplo concreto — **Testcase 2019 (ciclo=2019111500)**

| Tipo                 | Diretório base                                                       | Exemplos                                                                                      |
| -------------------- | -------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| **Análises**         | `${SUBMIT_HOME}/SMG/datainout/gsi/dataout/2019111500`                | `GANLCPT2019111500S.unf.TQ0299L064`                                                           |
| **Restarts**         | `${SUBMIT_HOME}/SMG/datainout/bam/pos/dataout/TQ0299L064/2019111500` | `GFCTCPT2019111500F.unf.TQ0299L064.convclP01`, `...outattP01`, `...outmdtP01`, `...sibprgP01` |
| **First Guess**      | `${SUBMIT_HOME}/SMG/datainout/bam/pos/dataout/TQ0299L064/2019111500` | `GANLCPT2019111500F.dir.TQ0299L064`, `GANLCPT2019111500F.fct.TQ0299L064`                      |
| **Previsões (GRIB)** | `${SUBMIT_HOME}/SMG/datainout/bam/pos/dataout/TQ0299L064/2019111500` | `GPOSCPT2019111500P.icn.TQ0299L064.grb`, `...inz...`, `GPOSCPT2019111509P.fct.TQ0299L064.grb` |
| **Prévia (GrADS)**   | `${SUBMIT_HOME}/SMG/datainout/bam/pre/dataout/`                      | `GANLNMC2019111500S.unf.TQ0299L064.GrADS`                                                     |

---

## 8. Dicas e solução de problemas

* **SVN “tag não encontrada”** → `svn list https://svn.cptec.inpe.br/smna/tag/` e escolha a mais recente.
* **Executável ausente** → ver `compile_smg.log` (palavras-chave: `error`, `fatal`, `not found`), conferir módulos/variáveis NETCDF/HDF5/PIO/PNETCDF.
* **Pós-processamento (GRIB) ausente** → habilitar no `run_cycle.sh` (trocar `"No"` → `"Yes"`).
* **Previsões 6h vs 3h** → `kpds=11` (6 h) ou `kpds=10` (3 h) no `POSTIN-GRIB`.
* **Rodar em background**:

  ```bash
  nohup ./run_cycle.sh ... > run_cycle.out &
  tail -f run_cycle.out
  ```

---

## Anexo: Início rápido (cheatsheet)

```bash
# 0) Obter o branch
ssh usuario@egeon.cptec.inpe.br
cd /home/${USER}
svn co https://svn.cptec.inpe.br/smna/branch/SMNA_v3.0.0.t11889
cd SMNA_v3.0.0.t11889/SMG

# 1) Ajustar paths
nano etc/mach/egeon_paths.conf
# nome_smg  SMNA_v3.0.0.t11889/SMG

# 2) Configurar (estrutura + fixos iniciais)
./config_smg.ksh configure -v -y

# 3) Compilar (sem ANGUPDATE)
./config_smg.ksh compile -v --all --no-ang

# 4) Verificar executáveis
./config_smg.ksh verify_executables -v --all --no-ang

# 5) Testcase para data atual (2025)
./config_smg.ksh copy_ncep_input 2025091500 -v

# 6) Pré-processamento (uma vez por resolução)
cd cptec/bam/run
./runPre -v -t 299 -l 64 -I 2025091500 -n 0 -p SMT -s -O -T -G

# 7) First Guess (9 h)
./runModel -das -np 128 -t 299 -l 64 -I 2025091500 -W 2025091509 -F 2025091509 -ts 3 -r -tr 6 -i 2 -p SMT -s sstwkl

# 8) Ciclo completo (exemplo)
cd ../../../run
./run_cycle.sh -t 299 -l 64 -gt 254 -p CPT -I 2025091506 -F 2025091618
```

---
