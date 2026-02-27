# Implementação Pré-Operacional do SMNA

> **Resumo:** Este guia explica, de ponta a ponta, como **obter**, **configurar**, **compilar** e **executar** o ciclo pré-operacional do **SMNA** (BAM + GSI), com foco inicial no **XC50** e preparação paralela para o **Egeon**.

## Sumário

* [1. Visão geral](#1-visão-geral)
* [2. Escopo e estratégia (XC50 ↔ Egeon)](#2-escopo-e-estratégia-xc50--egeon)
* [3. Obtenção do sistema (SVN)](#3-obtenção-do-sistema-svn)
* [4. Instalação e configuração](#4-instalação-e-configuração)
* [5. Compilação (GNU no XC50)](#5-compilação-gnu-no-xc50)
* [6. Execução: testcase e ciclo](#6-execução-testcase-e-ciclo)
* [7. Resultados e saídas geradas](#7-resultados-e-saídas-geradas)
* [8. Rodadas no Egeon (Tarefa #11598)](#8-rodadas-no-egeon-tarefa-11598)
* [9. Namelist de exemplo (MODELIN)](#9-namelist-de-exemplo-modelin)
* [10. Dicas e solução de problemas](#10-dicas-e-solução-de-problemas)

---

## 1. Visão geral

O **SMNA** (Sistema de Modelagem Numérica e Assimilação) integra:

* **BAM** (*Brazilian Atmospheric Model*): previsão numérica.
* **GSI** (*Gridpoint Statistical Interpolation*): assimilação de dados.

Para o ciclo se estabelecer corretamente, são necessários três “ingredientes”:

1. **Fluxo de dados**
2. **Estimativa *a priori*** do estado atmosférico (first guess)
3. **Estatísticas de erro** (observação e modelo)

Com base nisso, o plano pré-operacional adota as etapas:

* **Fluxo de dados**

  * Fase 1: uso de dados pré-processados pelo **NCEP**
  * Fase 2: uso de dados do **GTS** pré-processados pelo **INPE**
* **Estimativas *a priori***

  * Pré-operacionalização do **BAM** com janela de **9 h**
* **Estatísticas de erro**

  * Pré-operacionalização do **GSI**

---

## 2. Escopo e estratégia (XC50 ↔ Egeon)

Há duas possibilidades para a execução:

* **XC50**: alvo principal inicial (apesar de não ter garantia de hardware no momento);
* **Egeon**: receberá, em paralelo, os ajustes (scripts de submissão e build).

> **Decisão:** realizar o **pré-operacional completo no XC50**, enquanto os ajustes para o **Egeon** são implementados em paralelo.

---

## 3. Obtenção do sistema (SVN)

O código está no **SVN do CPTEC**, em dois sabores:

* **`trunk`**: desenvolvimento (para contribuidores)
* **`tag`**: versões estáveis (para usuários)

### Passo a passo (XC50)

1. Login no XC50:

```bash
ssh usuario@login-xc50.cptec.inpe.br -XC
```

2. Vá ao seu diretório no Lustre:

```bash
cd /lustre_XC50/${USER}
```

3. **Export** de uma versão publicada (exemplo `SMNA_v2.3.1`):

```bash
svn export https://svn.cptec.inpe.br/smna/tag/SMNA_v2.3.1
```

4. Para listar as versões/tags disponíveis:

```bash
svn list https://svn.cptec.inpe.br/smna/tag/
```

5. Para contribuir (controle de versão ativo), use **checkout** do `trunk`:

```bash
svn co https://svn.cptec.inpe.br/smna/trunk/SMNA
```

> **Resultado esperado:** um diretório como `/lustre_XC50/${USER}/SMNA_v2.3.1` (tags) ou `/lustre_XC50/${USER}/SMNA` (trunk) contendo todo o pacote.

---

## 4. Instalação e configuração

Entre na pasta do **SMG** (gerenciador do SMNA) e rode o *config*:

```bash
cd /lustre_XC50/${USER}/SMNA_v2.3.1/SMG
./config_smg.ksh
```

### Comandos do `config_smg.ksh`

* `configurar` – cria **estrutura de diretórios/links**
* `compilar` – **compila** GSI+BAM e utilitários
* `testcase` – executa um **caso de teste**
* `ajuda` – mostra **help** (executado também sem argumentos)

### Arquivo de paths (`etc/paths.conf`)

Revise os valores principais:

| Variável      | Valor sugerido (XC50)              |
| ------------- | ---------------------------------- |
| `nome_smg`    | `SMG`                              |
| `HOME`        | `/lustre_XC50/${USER}/SMNA_v2.3.1` |
| `SUBMIT_HOME` | `${HOME}`                          |

> **Atenção**
>
> * O pacote nasceu para a Tupã e usa `SUBMIT_HOME`/`WORK_HOME`. No **XC50**, aponte ambos para o **mesmo lugar** de `HOME`.
> * Se você alterar variáveis em `paths.conf`, **ajuste** os comandos deste guia conforme o seu ambiente.

### Criar a estrutura de diretórios

```bash
./config_smg.ksh configurar
```

Por padrão, a estrutura de **entrada/saída** e submissão é organizada sob `${SUBMIT_HOME}/SMNA` e `${WORK_HOME}/SMNA/SMG`. Se precisar personalizar nomes/locais, edite a função `vars_export` dentro de `config_smg.ksh` antes de executar `configurar`.

---

## 5. Compilação (GNU no XC50)

> **Recomendado:** compilar com o **toolchain GNU** no XC50.

### Carregar módulos

```bash
module load pbs
module load craype-x86-skylake
module load craype-network-aries
module load cray-netcdf
module swap PrgEnv-cray PrgEnv-gnu
```

### Compilar

1. Garanta o diretório correto:

```bash
cd /lustre_XC50/${USER}/SMNA_v2.3.1/SMG
```

2. Configure ambiente e estrutura (se ainda não fez):

```bash
./config_smg.ksh configurar
```

3. Compile e salve o log:

```bash
./config_smg.ksh compilar 2>&1 | tee compile_smg.log
```

4. Acompanhe em tempo real:

```bash
tail -f compile_smg.log
```

### Verificar executáveis gerados

```bash
# GSI
ls /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/cptec/bin/gsi.exe
ls /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/cptec/bin/gsi_angupdate.exe

# IncTime
ls /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/cptec/bin/inctime

# Pré-processamento (ver lista esperada na wiki interna)
ls /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/cptec/bam/pre/exec/

# BAM
ls /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/cptec/bam/model/exec/ParModel_MPI

# Pós-processamento
ls /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/cptec/bam/pos/exec/PostGrib
```

> Se algo **não existir**, abra o `compile_smg.log`, **identifique e corrija** o erro, e repita a compilação.

---

## 6. Execução: testcase e ciclo

O **testcase** executa GSI+BAM para **TQ0062L028 (\~200 km)** e **TQ0254L064 (\~50 km)** gerando **análises** e **previsões de 9 h** (saídas horárias) para datas de **jan/2013** ou **mai/2015**.
Os **dados** ficam em:
`/scratchin/grupos/assim_dados/home/gdad/public` (leia o `README` lá).

> **Obs.:** por padrão, esta versão assimila **dados convencionais** (PrepBUFR/BUFR).

### 6.1 Preparação

1. No diretório do SMG:

```bash
cd /lustre_XC50/${USER}/SMNA_v2.3.1/SMG
./config_smg.ksh copy_fixed_files
```

2. Linkar dados necessários (exemplo):

```bash
cd /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/datainout/bam/pre/datain
ln -s /lustre_xc50/ioper/data/external/2021010100/dataout/Umid_Solo/GL_SM.GPNR.2021010100.vfm
ln -s /lustre_xc50/ioper/data/external/2021010100/dataout/NCEP/rtgssthr_grb_0.083.grib2.20210101
ln -s /lustre_xc50/ioper/data/external/2021010100/dataout/NCEP/gblav.T00Z.atmanl.nemsio.2021010100
```

### 6.2 Pré-processamento do BAM (uma única vez)

```bash
cd /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/cptec/bam/run
/bin/bash runPre -v -t 299 -l 64 -I 2021010100 -n 0 -p SMT -s -O -T -G -Gp gblav -Gt Grid
# sem cores:
# /bin/bash runPre ... | sed $'s/\e\[[0-9;:]*[a-zA-Z]//g'
```

### 6.3 First Guess (9 h) e restarts

```bash
/bin/bash ./runModel -das -np 480 -N 10 -d 4 -t 299 -l 64 -I 2021010100 -W 2021010109 -F 2021010109 -ts 3 -r -tr 6 -i 2 -p SMT -s sstwkl
# sem cores:
# /bin/bash ./runModel ... | sed $'s/\e\[[0-9;:]*[a-zA-Z]//g'
```

### 6.4 Ciclo de assimilação

1. Ajuda/uso:

```bash
cd /lustre_XC50/${USER}/SMNA_v2.3.1/SMG/run
./run_cycle.sh -h
```

2. Execução típica:

```bash
./run_cycle.sh -t 299 -l 64 -gt 254 -p CPT -I 2021010106 -F 2021010818
```

3. Execução com `nohup` (reconectável) e log limpo:

```bash
nohup ./run_cycle.sh -t 299 -l 64 -gt 254 -p CPT -I 2021010106 -F 2021010818 \
  | sed $'s/\e\[[0-9;:]*[a-zA-Z]//g' > run_cycle.out &
tail -f run_cycle.out
```

> **Notas úteis**
>
> * **`-np 480`** no BAM: número de processos coincide com a **partição de restarts** (superfície/atmosfera/radiação/convecção/nuvens).
> * O **pós-processamento** do BAM vem **desativado** no ciclo. Para habilitar, no `run_cycle.sh` mude a chamada final do `run_model.sh` de `"No"` para `"Yes"` (gera GRIB).
> * Para previsões a **cada 6 h**, mude `kpds=10` para `kpds=11` em `POSTIN-GRIB` (ver pasta do `run`).

---

## 7. Resultados e saídas geradas

### Análises

* **Dir:** `/lustre_XC50/${USER}/SMNA_v2.3.1/SMG/datainout/gsi/dataout/<DATA>`
* **Arquivo:** `GANLCPT<AAAAMMDD><HH>S.unf.TQ0062L028`
* **Uso:** arquivo de análise lido pelo BAM (com estatísticas e logs do GSI no mesmo diretório).

### Restarts

* **Dir:** `/lustre_XC50/${USER}/SMNA_v2.3.1/SMG/datainout/bam/pos/dataout/<RESOLUCAO>/<DATA>`
* **Arquivos típicos:**

  * `GFCTCPT<ana><ana>F.unf.<res>.convclP<rank>` (convecção/radiação/nuvens)
  * `GFCTCPT<ana><ana>F.unf.<res>.outattP<rank>` (dinâmica – tempo corrente)
  * `GFCTCPT<ana><ana>F.unf.<res>.outmdtP<rank>` (dinâmica – tempo passado)
  * `GFCTCPT<ana><ana>F.unf.<res>.sibprgP<rank>` (superfície)

### First Guess

* **Dir:** `/lustre_XC50/${USER}/SMNA_v2.3.1/SMG/datainout/bam/pos/dataout/<RESOLUCAO>/<DATA>`
* **Arquivos:**

  * `GANLCPT<ana><prev>F.dir.<res_mcga>` (header ASCII)
  * `GANLCPT<ana><prev>F.fct.<res_mcga>` (espectral)

### Previsões (GRIB)

* **Dir:** `/lustre_XC50/${USER}/SMNA_v2.3.1/SMG/datainout/bam/pos/dataout/<RESOLUCAO>/<DATA>`
* **Arquivos:**

  * `GPOSCPT<ana><ana>P.icn.<res_mcga>.grb` (análise pós-processada)
  * `GPOSCPT<ana><ana>P.inz.<res_mcga>.grb` (análise inicializada – modos normais)
  * `GPOSCPT<ana><prev>P.fct.<res_mcga>.grb` (previsão pós-processada)

> **Prévia em grade (full):** ver `GANLNMC<ana>S.unf.<res_mcga>.GrADS` em
> `/scratchout/grupos/<grupo>/home/<usuario>/SMG/datainout/bam/pre/dataout`.

---

## 8. Rodadas no Egeon (Tarefa #11598)

### Ambiente (Intel + MPI) utilizado

```bash
#!/bin/bash -x
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

### Recursos SLURM (exemplo)

```bash
#SBATCH --ntasks=160
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
```

* No **T299L064**, \~**15 dias** de previsão em \~**30 min** (referência empírica do ambiente).
* Operacional: **T666L064**, \~**11 dias** em \~**140 min**.

> **Observação de custo:** o aumento de custo entre T299→T666 é consistente com crescimento \~quadrático, dominado pelas transformadas (Fourier/Legendre).

Exemplos de fila:

```
174222  batch  BAM_1003  ioper  R  2:20:35  10 n[14-23]
154466  batch  bam enver.ra  R  31:23  10 n[03-04,06-07,09-13,28]
```

---

## 9. Namelist de exemplo (MODELIN)

<details><summary>Clique para expandir o MODELIN (exemplo T299L064)</summary>

```fortran
&MODEL_RES
 trunc    =0299,
 vert     =064,
 dt       =450.0,
 ...
/
&MODEL_IN
 slagr    =.TRUE.,
 ...
/
&PHYSPROC
 UNIFIED  = .TRUE.
 ISWRAD   = 'CRD'
 ILWRAD   = 'CRD'
 ICCON    = 'ARA'
 ISCON    = 'TIED'
 ILCON    = 'HUMO'
 IGWD     = 'YES'
 CRDCLD   = 6,
 ISIMP    = 'NO ',
 thermcell=1,
 atmpbl   = 4,
 schemes  = 3,
 ...
/
&PHYSCS
 swint   = 3600.0
 trint   = 10800.0
 co2val  = 370.0,
/
&COMCON
 initlz = 2,
 fint   = 6,
 ifsst  = -1,
 ...
/
&NumberOutPutForecast
  maxtfm2= 17
  cthrfx( 1: 4) =   6.0, 12.0, 18.0, 24.0,
  ...
/
```

</details>

> **Dica:** mantenha um **MODELIN-control** versionado e derive variações por `include`/patch (evita drift entre casos).

---

## 10. Dicas e solução de problemas

* **SVN “tag não encontrada”**

  * Rode `svn list https://svn.cptec.inpe.br/smna/tag/` e escolha a **última versão**.
* **Executável ausente após build**

  * Verifique `compile_smg.log` (busque por `error`, `fatal`, `not found`).
  * Confirme **módulos** e **variáveis de ambiente** de bibliotecas (NetCDF/HDF5/PIO/PNETCDF).
* **Diferença de caminhos (`/lustre_XC50` vs `/lustre_xc50`)**

  * No XC50, use **`/lustre_XC50`** (com `X` maiúsculo) para manter consistência.
* **Pós-processamento (GRIB) ausente**

  * Habilite na chamada do `run_model.sh` **(“Yes”)** dentro do `run_cycle.sh`.
* **Previsões 6h vs 3h**

  * Ajuste `kpds=11` (6h) ou `kpds=10` (3h) em `POSTIN-GRIB`.
* **Rodando com `nohup`**

  * Acompanhe com `tail -f run_cycle.out`; o processo **segue** mesmo após logout.

---

**Contato/Manutenção:** Atualize este documento conforme mudanças em scripts, caminhos, módulos e políticas de execução (XC50/Egeon).
