#!/usr/bin/env bash
#-----------------------------------------------------------------------------#
# Git Release Flow Script                                                     #
#-----------------------------------------------------------------------------#
# DESCRIPTION
#   Este script automatiza um fluxo simples e seguro de release em repositórios
#   Git. Ele cria uma tag de versão a partir de uma branch específica e depois
#   integra essa branch na branch principal do projeto.
#
#   O fluxo executado é o seguinte:
#
#       1) Checkout da branch de origem da release
#       2) Atualização da branch local com o repositório remoto
#       3) Criação de uma tag anotada da release
#       4) Push da tag para o repositório remoto
#       5) Merge da branch de origem na branch principal (main)
#       6) Push da branch principal atualizada
#       7) (Opcional) criação de uma nova branch de desenvolvimento
#
#   Esse processo ajuda a manter um histórico claro de releases e evita erros
#   manuais comuns ao executar esses passos separadamente.
#
# USAGE
#
#   git_release_flow.sh \
#       --source-branch <branch_origem> \
#       --tag <tag_release> \
#       [--next-branch <nova_branch_dev>]
#
# EXAMPLES
#
#   Criar release a partir da branch "release-v1":
#
#       ./git_release_flow.sh \
#           --source-branch release-v1 \
#           --tag v1.0.0
#
#   Criar release e iniciar nova branch de desenvolvimento:
#
#       ./git_release_flow.sh \
#           --source-branch release-v1 \
#           --tag v1.0.0 \
#           --next-branch develop-v1.1
#
# WHEN TO USE
#
#   Este script deve ser utilizado quando:
#
#   - Uma versão do software está pronta para ser publicada.
#   - Deseja-se marcar oficialmente essa versão com uma tag Git.
#   - A branch de desenvolvimento precisa ser integrada ao branch principal.
#
# WHY THIS SCRIPT EXISTS
#
#   Criar releases manualmente envolve vários comandos Git que podem gerar
#   inconsistências se executados incorretamente. Este script garante que
#   o fluxo seja executado sempre na mesma ordem e com confirmação explícita
#   do usuário antes de modificar o repositório.
#
# SAFETY FEATURES
#
#   - set -euo pipefail evita execução silenciosa de erros
#   - confirmação interativa antes de executar operações destrutivas
#   - uso de merge --no-ff para preservar histórico de branches
#
# REQUIREMENTS
#
#   - git instalado
#   - repositório já clonado
#   - acesso de escrita ao remoto (origin)
#
#-----------------------------------------------------------------------------#

set -euo pipefail

#-----------------------------------------------------------------------------
# Default configuration
#-----------------------------------------------------------------------------
# Branch principal do projeto. Pode ser alterada aqui caso o repositório
# utilize outro nome como "master" ou "stable".
MAIN_BRANCH="main"

#-----------------------------------------------------------------------------
# Argument variables
#-----------------------------------------------------------------------------
# SOURCE_BRANCH : branch de onde a release será criada
# TAG           : tag da release (ex: v1.2.0)
# NEXT_BRANCH   : branch opcional para iniciar o próximo ciclo de desenvolvimento
SOURCE_BRANCH=""
TAG=""
NEXT_BRANCH=""

#-----------------------------------------------------------------------------
# Argument parsing
#-----------------------------------------------------------------------------
# Processa os argumentos passados na linha de comando.
while [[ $# -gt 0 ]]; do
    case "$1" in
        --source-branch)
            SOURCE_BRANCH="$2"
            shift 2
            ;;
        --tag)
            TAG="$2"
            shift 2
            ;;
        --next-branch)
            NEXT_BRANCH="$2"
            shift 2
            ;;
        *)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done

#-----------------------------------------------------------------------------
# Argument validation
#-----------------------------------------------------------------------------
# SOURCE_BRANCH e TAG são obrigatórios para o fluxo de release.
[[ -n "$SOURCE_BRANCH" && -n "$TAG" ]] || {
    echo "Usage:"
    echo "  git_release_flow.sh --source-branch BRANCH --tag TAG [--next-branch BRANCH]"
    exit 1
}

#-----------------------------------------------------------------------------
# Show configuration summary
#-----------------------------------------------------------------------------
# Mostra ao usuário o que será feito antes de executar.
echo "[INFO] source branch : $SOURCE_BRANCH"
echo "[INFO] tag           : $TAG"
echo "[INFO] next branch   : $NEXT_BRANCH"
echo

#-----------------------------------------------------------------------------
# Confirmation prompt
#-----------------------------------------------------------------------------
# Solicita confirmação explícita antes de executar comandos Git que alteram
# o histórico do repositório.
read -p "Proceed? [y/N]: " ans
[[ "$ans" =~ ^[Yy]$ ]] || exit 0

#-----------------------------------------------------------------------------
# Checkout source branch
#-----------------------------------------------------------------------------
# Garante que estamos na branch correta para gerar a release.
echo "[INFO] checkout source branch"
git checkout "$SOURCE_BRANCH"

#-----------------------------------------------------------------------------
# Update branch
#-----------------------------------------------------------------------------
# Sincroniza com o remoto para evitar que a tag seja criada sobre uma
# versão desatualizada da branch.
echo "[INFO] pull latest"
git pull

#-----------------------------------------------------------------------------
# Create release tag
#-----------------------------------------------------------------------------
# Cria uma tag anotada, que inclui mensagem e metadata.
echo "[INFO] creating tag"
git tag -a "$TAG" -m "Release $TAG"

#-----------------------------------------------------------------------------
# Push tag
#-----------------------------------------------------------------------------
# Publica a tag no repositório remoto.
echo "[INFO] pushing tag"
git push origin "$TAG"

#-----------------------------------------------------------------------------
# Merge into main
#-----------------------------------------------------------------------------
# Integra a branch da release na branch principal do projeto.
# O uso de --no-ff garante que o histórico preserve a existência da branch.
echo "[INFO] merging into $MAIN_BRANCH"
git checkout "$MAIN_BRANCH"
git pull
git merge --no-ff "$SOURCE_BRANCH"

#-----------------------------------------------------------------------------
# Push updated main branch
#-----------------------------------------------------------------------------
echo "[INFO] pushing main"
git push origin "$MAIN_BRANCH"

#-----------------------------------------------------------------------------
# Optional next development branch
#-----------------------------------------------------------------------------
# Caso especificado, cria automaticamente a próxima branch de desenvolvimento.
if [[ -n "$NEXT_BRANCH" ]]; then
    echo "[INFO] creating next development branch $NEXT_BRANCH"
    git checkout -b "$NEXT_BRANCH"
    git push origin "$NEXT_BRANCH"
fi

#-----------------------------------------------------------------------------
# Completion message
#-----------------------------------------------------------------------------
echo "[OK] release flow completed"
