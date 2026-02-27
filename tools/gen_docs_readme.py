#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from datetime import datetime

REPO_ROOT = Path(__file__).resolve().parents[1] if (Path(__file__).name != "<stdin>") else Path.cwd()
DOCS_DIR = REPO_ROOT / "docs"
OUT = DOCS_DIR / "README.md"

def is_hidden(p: Path) -> bool:
    return any(part.startswith(".") for part in p.parts)

def md_link(path: Path) -> str:
    # links relativos ao docs/
    rel = path.relative_to(DOCS_DIR).as_posix()
    return f"- [`{rel}`]({rel})"

def main() -> int:
    if not DOCS_DIR.exists() or not DOCS_DIR.is_dir():
        raise SystemExit(f"[ERRO] Diretório não encontrado: {DOCS_DIR}")

    # coleta (1 nível de profundidade: docs/*)
    entries = sorted([p for p in DOCS_DIR.iterdir() if not is_hidden(p) and p.name != "README.md"],
                     key=lambda p: (0 if p.is_dir() else 1, p.name.lower()))

    dirs = [p for p in entries if p.is_dir()]
    files = [p for p in entries if p.is_file()]

    lines: list[str] = []
    lines.append("# Documentação (docs/)\n")
    lines.append("Este diretório reúne documentação técnica do **SMNA** (BAM + GSI) e materiais auxiliares.\n")
    lines.append("> Este índice é gerado automaticamente. Se tu adicionar/remover arquivos em `docs/`, rode o gerador novamente.\n")
    lines.append(f"_Última atualização do índice: {datetime.now().strftime('%Y-%m-%d %H:%M')}_\n")

    if dirs:
        lines.append("## Pastas\n")
        for d in dirs:
            lines.append(md_link(d))
        lines.append("")  # linha em branco

    if files:
        lines.append("## Arquivos\n")
        for f in files:
            lines.append(md_link(f))
        lines.append("")

    # dica opcional
    lines.append("## Como atualizar este índice\n")
    lines.append("```bash\npython3 tools/gen_docs_readme.py\n```\n")

    OUT.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"[OK] Gerado: {OUT}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
