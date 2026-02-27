# 🗂️ Diretórios Equivalentes entre Clusters

Este documento lista a correspondência entre os principais diretórios de trabalho utilizados nos clusters **BASTOS**, **EGEON** e **XC50**.

| BASTOS              | EGEON                 | XC50 |
|---------------------|-----------------------|------|
| `/dados/das`        | `/pesq/dados/das`     | —    |
| `/share/das`        | `/pesq/share/das`     | —    |
| —   | `/mnt/beegfs/home/joao.gerd` | `/lustre_xc50/joao_gerd` |


---

### 🧭 Observações

- O cluster **XC50** utiliza o sistema de arquivos **Lustre**, com prefixo `/lustre_xc50/`.
- O cluster **EGEON** utiliza o sistema de arquivos **beegfs**, com prefixo `/mnt/beegfs`.  
- Caminhos podem variar conforme o grupo ou usuário.  
- Recomenda-se manter scripts compatíveis usando variáveis como:

```bash
  export BASE_DIR=${BASEDIR:-/scratchin/grupos/assim_dados}
````

---

📄 **Autor:** João Gerd Zell de Mattos
📅 **Atualizado:** 13/10/2025
📘 **Licença:** CC BY-NC 3.0



