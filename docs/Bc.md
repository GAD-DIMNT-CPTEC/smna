# Guia prático: Spin-up de coeficientes de correção de viés (satbias) no GSI e seleção de canais

> Este guia resume e operacionaliza as recomendações do manual de Assimilação de Radiâncias do GSI para casos em que você ainda não possui coeficientes de viés “maduros” (mass e angle-dependent) e vai rodar janelas curtas (ex.: ~1 semana). Inclui: pré-requisitos, checagens de consistência (satinfo × satbias), procedimento de spin-up e um script Bash de exemplo para automatizar o processo.




---

1) Visão geral

O GSI distribui arquivos de exemplo em ./fix/:

satbias_in: sample.satbias

satbias_angle: global_satangbias.txt e nam_global_satangbias.txt


> Atenção: são amostras. Cada usuário deve gerar seus próprios coeficientes, e normalmente é preciso semanas/meses de ciclo para estabilizá-los.



Para experimentos curtos, a orientação prática é “antecipar” um spin-up rápido:

1. Iniciar com coeficientes de uma data o mais próxima possível do seu caso.


2. Rodar uma análise do GSI com mass bias e angle-dependent habilitados → isso atualiza os arquivos de coeficientes.


3. Repetir a mesma análise, com o mesmo background e observações, porém alimentando os coeficientes atualizados.


4. Repetir o passo 3 ~10 vezes (spin-up dos coeficientes).


5. Passar para o próximo ciclo/data e repetir (passos 2–4).


6. Após 1–2 dias, o mass bias tende a estabilizar; angle-dependent converge mais lentamente.





---

2) Pré-requisitos e boas práticas

Arquivos chave

satinfo (lista de canais e usage flags)

satbias_in (coeficientes de viés “mass” por canal/sensor)

satbias_angle (coeficientes dependentes de ângulo)


Consistência de canais

Canais em satinfo devem bater com os de satbias_in e satbias_angle.

Se um canal existe em satinfo e não em satbias_in, o GSI adiciona no update com 0 como valor inicial.

Se um canal existe em satbias_in e não em satinfo, ele é removido do arquivo atualizado.

Se satinfo × satbias_angle não combinam, o GSI usa o que tiver para a análise, mas a ferramenta de atualização de ângulo pode falhar. Garanta a coincidência antes de rodar a atualização de ângulo.


Seleção de canais (além de bias)

1. Topo do modelo × função peso do instrumento: se a função peso de um canal ultrapassa o topo do seu modelo, exclua esse canal.


2. Canais problemáticos (calibração/defeito): desligue (monitor apenas).


3. Testes de impacto: monitore/avalie por um período e promova de “monitor” → “use” apenas após verificar impacto positivo.





---

3) Fluxo recomendado de spin-up (execução curta)

Para cada janela de análise (ex.: YYYYMMDDHH):

1. Preparar coeficientes iniciais

Se não tem coeficientes próprios: copie os mais próximos no tempo (ou use os samples como último recurso).



2. Análise 1 (update inicial)

Rode o GSI com mass e angle habilitados → gera satbias_out e arquivos de atualização de ângulo.



3. Iterações de spin-up (≈ 10 vezes)

Para cada iteração: use o mesmo background e observações da Análise 1, porém alimente os coeficientes atualizados da iteração anterior.



4. Avançar para a próxima data

Repita passos 2–3.



5. Após 1–2 dias

Use o mass bias spin-ado em seus casos reais. Continue girando angle-dependent se necessário.





---

4) Checagens automáticas de consistência de canais

Antes de rodar, valide:

satinfo vs satbias_in → para evitar “canais órfãos” ou adicionados como zero indevidamente.

satinfo vs satbias_angle → para não quebrar a ferramenta de angle update.


O script abaixo inclui uma função que lista e compara os conjuntos de canais.


---

5) Script Bash de exemplo para spin-up rápido

> Objetivo: automatizar o spin-up de satbias_in e satbias_angle por data de análise, com N iterações usando mesmo BG/OBS, atualizando os coeficientes a cada loop.
Pré-condições: você já tem um wrapper runGSI que:

lê satinfo, satbias_in, satbias_angle do diretório de trabalho;

produz satbias_out e arquivos de ângulo atualizados (nomenclatura pode variar no seu sistema).


```bash

#!/usr/bin/env bash
set -euo pipefail

# ===================== Configuração do usuário ======================
# Datas de análise (ex.: 2 dias + período alvo)
AN_DATES=("20250128_12" "20250128_18" "20250129_00" "20250129_06")

# Número de iterações de spin-up por data (recomenda-se ~10)
SPINUP_ITERS=10

# Caminhos de entrada
FIX_DIR="/path/para/fix"                         # contém samples e/ou coeficientes base
SATINFO="/path/para/fix/satinfo"                 # satinfo mestre
SATBIAS_SEED="${FIX_DIR}/sample.satbias"         # ponto de partida (ideal: mais próximo da data)
SATANG_SEED="${FIX_DIR}/global_satangbias.txt"   # idem (ou nam_global_satangbias.txt)

# Executável / wrapper do GSI
RUN_GSI="/caminho/para/runGSI"                   # seu script que roda o GSI nesta data

# Diretório raiz de trabalho (um subdir por data)
WORK_ROOT="/path/work/gsi_spinup"

# ===================== Funções auxiliares ===========================

err() { echo "[ERROR] $*" >&2; exit 1; }
ok()  { echo "[OK] $*"; }
inf() { echo "[INFO] $*"; }

# Extrai lista de "canais" (definição simplificada) dos arquivos
# Ajuste o awk conforme o formato real da sua instalação.
list_channels_satinfo() {
  awk 'NF>0 && $1 !~ /^#/ {print $1}' "$1" | sort -u
}
list_channels_satbias_in() {
  awk 'NF>0 && $1 !~ /^#/ {print $1}' "$1" | sort -u
}
list_channels_satang() {
  awk 'NF>0 && $1 !~ /^#/ {print $1}' "$1" | sort -u
}

# Compara conjuntos e avisa diferenças
check_consistency() {
  local satinfo="$1" satbias="$2" satang="$3"

  inf "Checando consistência de canais (satinfo × satbias_in × satbias_angle)..."
  comm -3 <(list_channels_satinfo "$satinfo") <(list_channels_satbias_in "$satbias") \
    | sed 's/^/\t/;' || true

  # Verificação satinfo × satang
  local tmp_sf tmp_sa
  tmp_sf=$(mktemp); tmp_sa=$(mktemp)
  list_channels_satinfo "$satinfo" > "$tmp_sf"
  list_channels_satang  "$satang"  > "$tmp_sa"

  if ! diff -u "$tmp_sf" "$tmp_sa" >/dev/null 2>&1; then
    echo "[WARNING] satinfo e satbias_angle NÃO coincidem 1-a-1."
    echo "          O GSI pode rodar, mas o ATUALIZADOR DE ÂNGULO pode falhar."
    echo "          Revise os canais antes de atualizar ângulo."
  else
    ok "satinfo e satbias_angle compatíveis."
  fi
  rm -f "$tmp_sf" "$tmp_sa"
}

# Copia com progresso quando rsync existir
copy_file() {
  local src="$1" dst="$2"
  mkdir -p -- "$(dirname -- "$dst")"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --human-readable --info=progress2,name0 -i -- "$src" "$dst"
  else
    cp -pf -- "$src" "$dst"
  fi
}

# Roda uma análise do GSI em um diretório de trabalho já populado
run_one_analysis() {
  local wdir="$1"
  ( cd "$wdir" && "$RUN_GSI" )
}

# Após a análise, promove satbias_out → satbias_in e atualiza ângulo
promote_coeffs() {
  local wdir="$1"
  local out="${wdir}/satbias_out"
  local in="${wdir}/satbias_in"
  local ang_in="${wdir}/satbias_angle"
  local ang_upd="${wdir}/satbias_angle.updated"   # ajuste ao nome que seu fluxo gera

  [[ -s "$out" ]] || err "satbias_out não encontrado/ vazio em $wdir"
  copy_file "$out" "$in"
  ok "Atualizado: satbias_in <= satbias_out"

  if [[ -s "$ang_upd" ]]; then
    copy_file "$ang_upd" "$ang_in"
    ok "Atualizado: satbias_angle <= satbias_angle.updated"
  else
    echo "[NOTICE] Arquivo de ângulo atualizado não encontrado; mantendo o atual."
  fi
}

# ===================== Pipeline principal ===========================

mkdir -p -- "$WORK_ROOT"

for AD in "${AN_DATES[@]}"; do
  WDIR="${WORK_ROOT}/${AD}"
  inf "Preparando diretório de trabalho: $WDIR"
  mkdir -p -- "$WDIR"

  # 1) Seed inicial dos coeficientes (mais perto no tempo do caso real)
  copy_file "$SATBIAS_SEED" "${WDIR}/satbias_in"
  copy_file "$SATANG_SEED"  "${WDIR}/satbias_angle"
  copy_file "$SATINFO"      "${WDIR}/satinfo"

  # 2) Checagem de consistência
  check_consistency "${WDIR}/satinfo" "${WDIR}/satbias_in" "${WDIR}/satbias_angle"

  # 3) Análise 1 (gera atualização inicial)
  inf "Rodando análise inicial (update de coeficientes)..."
  run_one_analysis "$WDIR"
  promote_coeffs "$WDIR"

  # 4) Iterações de spin-up (mesmo BG/OBS, só trocando coeficientes atualizados)
  for iter in $(seq 1 "$SPINUP_ITERS"); do
    inf "Spin-up ${iter}/${SPINUP_ITERS} — data ${AD}"
    run_one_analysis "$WDIR"
    promote_coeffs "$WDIR"
  done

  ok "Concluído spin-up para ${AD}"
done

ok "Pipeline finalizado. Considere usar os coeficientes 'mass' após 1–2 dias e continuar girando 'angle' se necessário."
```
Notas sobre o script

Idempotência: se você reexecutar para a mesma data, o diretório é reaproveitado e os coeficientes são promovidos a cada iteração.

Formato dos arquivos: os awk de extração de “canais” são genéricos; adeque ao layout real dos seus satinfo/satbias.

Atualização de ângulo: o nome do arquivo de saída pode variar no seu fluxo (ajuste satbias_angle.updated).

runGSI: o wrapper deve respeitar satbias_in, satbias_angle e satinfo do diretório de trabalho.



---

6) Dicas de operação

Escolha dos canais: antes de “ligar” um canal:

Cheque a função peso vs topo do modelo.

Monitore o canal por um período; acompanhe a série temporal do viés.

Promova de monitor → uso somente após avaliar o impacto na previsão.


Falhas na atualização de ângulo:

Quase sempre relacionadas a mismatch satinfo × satbias_angle.

Corrija a lista de canais e retente.


Janelas curtas:

Faça o spin-up 2 dias antes do período-alvo para chegar com mass razoavelmente estabilizado.




---

7) Estrutura de diretórios sugerida
```
project/
├─ fix/
│  ├─ satinfo
│  ├─ sample.satbias
│  ├─ global_satangbias.txt
│  └─ nam_global_satangbias.txt
├─ scripts/
│  ├─ runGSI              # seu wrapper
│  └─ gsi_spinup.sh       # (este script)
└─ work/
   ├─ 20250128_12/
   │  ├─ satinfo
   │  ├─ satbias_in
   │  ├─ satbias_out
   │  ├─ satbias_angle
   │  └─ satbias_angle.updated
   └─ 20250128_18/
      └─ ...

```
---

8) Checklist rápido

[ ] satinfo consistente com satbias_in e satbias_angle.

[ ] Spin-up com mesmo BG/OBS por ~10 iterações/ciclo.

[ ] Avançar ciclo e repetir.

[ ] Após 1–2 dias, usar mass bias nos casos reais; continuar girando angle.

[ ] Seleção de canais: função peso × topo do modelo; status de calibração; impacto.

