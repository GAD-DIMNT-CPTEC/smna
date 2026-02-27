# PostAlt — Passo-a-passo de Compilação (Local, Egeon e XC50)

Este README reúne **o passo-a-passo completo** para compilar o **PostAlt**, incluindo:
- Checkout das dependências (**sharedLibs**)
- Compilação de **libmisc** e **sigioBAM**
- **Autotools** do PostAlt (`autogen.sh` + `configure`)
- Diferenças de ambiente para **máquina local**, **Egeon (SLURM)** e **XC50 (PBS/Cray)**

---

## 1) Pré-requisitos

### 1.1 Ferramentas básicas
- `svn`, `make`, `autoconf`, `automake`, `libtool`, `pkg-config`
- Compiladores C e Fortran:
  - **Local**: `gcc`, `gfortran` (ou `clang`/`ifort` se preferir)
  - **Egeon**: módulos `gcc`/`openmpi` (ou `intel`)
  - **XC50**: wrappers `cc`/`ftn` via `PrgEnv-*`

### 1.2 Estrutura de diretórios sugerida
```bash
$HOME/
  PostAlt/               # (código-fonte do PostAlt)
  sharedLibs/            # (obtido via SVN)
    libmisc/
    libsigiobam/
````

---

## 2) Variáveis que ajudam o `configure` (opcional mas recomendável)

Use quando precisar **apontar explicitamente** onde estão `include/` e `lib/` do `sharedLibs` ou ajustar compiladores:

```bash
# Paths para headers e libs das deps:
export CPPFLAGS="-I$HOME/sharedLibs/libmisc/include -I$HOME/sharedLibs/libsigiobam/include"
export CFLAGS="$CFLAGS"
export FCFLAGS="$FCFLAGS"
export LDFLAGS="-L$HOME/sharedLibs/libmisc/lib -L$HOME/sharedLibs/libsigiobam/lib"

# Ordem de link:
export LIBS="-lsigiobam -lmisc"

# Escolher compiladores conforme ambiente:
# Local (exemplo):
export CC=gcc
export FC=gfortran

# Egeon (exemplo):
# export CC=gcc; export FC=gfortran

# XC50 (Cray):
# export CC=cc;  export FC=ftn
```

> 💡 Observação: a sintaxe tradicional do `configure` é `--prefix=$HOME/...`.

---

## 3) Ambiente por máquina (módulos)

### 3.1 **Egeon** (SLURM)

Carregar módulos (exemplo genérico):

```bash
module purge
module load gnu9/9.4.0
module load openmpi4/4.1.1

```

### 3.2 **XC50** (PBS / Cray)

```bash
module purge
module load PrgEnv-gnu          # ou PrgEnv-intel

```

---

## 4) Passo-a-passo oficial

### 4.1 Baixe o `sharedLibs` e compile `libmisc` e `sigioBAM`

```bash
cd $HOME
svn co https://svn.cptec.inpe.br/slib/trunk/sharedLibs

cd $HOME/sharedLibs/libmisc
./autogen.sh

# se estiver na tua máquina local:
./configure --prefix=$HOME/sharedLibs

# se estiver na egeon:
./configure --prefix=$HOME/sharedLibs --enable-egeon

# se estiver na XC50:
./configure --prefix=$HOME/sharedLibs --enable-cray

make
make install
```

```bash
cd $HOME/sharedLibs/libsigiobam

# se estiver na tua máquina local:
./configure --prefix=$HOME/sharedLibs

# se estiver na egeon:
./configure --prefix=$HOME/sharedLibs --enable-egeon

# se estiver na XC50:
./configure --prefix=$HOME/sharedLibs --enable-cray

make
make install
```

### 4.2 Agora baixe o PostAlt

```bash
cd $HOME
svn co https://svn.cptec.inpe.br/ad/trunk/PostAlt

cd postAlt
./autogen.sh

# se estiver na tua máquina local:
./configure --prefix=$HOME

# se estiver na egeon:
./configure --prefix=$HOME --enable-egeon

# se estiver na XC50:
./configure --prefix=$HOME --enable-cray

make
make install
```

> Deve ter sido criado um diretório **`$HOME/bin`**. Dentro dele deve estar o **`postAlt`**.

---

## 5) Teste rápido

```bash
$HOME/bin/postAlt
```

Executa **sem argumentos** para mostrar as instruções de uso.

---
## 6) Caso ocorra erro no `configure`

Se surgir algo como:

```
configure: error: miscellany library (libmisc) was not found, try specifying --with-misc
```

Execute novamente o `configure` informando explicitamente os caminhos das bibliotecas:

```bash
./configure --prefix=$HOME \
            --with-misc=$HOME/sharedLibs \
            --with-sigiobam=$HOME/sharedLibs
```

> Isso indica que o script não localizou automaticamente `libmisc` e `libsigiobam`.
> O comando acima aponta manualmente para os diretórios onde as bibliotecas foram instaladas.

---

## 7) Dicas rápidas

* **Linkagem**: se aparecer “undefined reference”, confira a **ordem das libs** em `LIBS` e os caminhos em `LDFLAGS`.
* **Toolchain**: compile tudo com o **mesmo compilador** (ex.: todos com `gcc/gfortran` ou `cc/ftn`).
* **Módulos**: em builds via SLURM/PBS, carregue os mesmos módulos no *job* que usou no *login*.

---

**Autor:** João Gerd Zell de Mattos
**Instituição:** CPTEC/INPE
**Licença:** CC BY-NC 3.0
