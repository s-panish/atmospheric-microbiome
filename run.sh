#!/usr/bin/env bash
# Точки запуска только две: amplicon_16S_ITS.Rmd (по одному маркеру за раз), shotgun.Rmd.
set -euo pipefail
cd "$(dirname "$0")"
mode="${1:-}"
config="${2:-config.yaml}"
if [[ ! -f "$config" ]]; then echo "Конфигурация не найдена: $config" >&2; exit 2; fi
render_amplicon() {
  local marker="$1"
  mkdir -p "results/$marker"
  Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("amplicon_16S_ITS.Rmd", params=list(config=a[1],marker=a[2]), output_file="report.html", output_dir=file.path("results",a[2]), envir=new.env(parent=globalenv()))' "$config" "$marker"
}
render_shotgun() {
  mkdir -p results/shotgun
  Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("shotgun.Rmd", params=list(config=a[1]), output_file="report.html", output_dir="results/shotgun", envir=new.env(parent=globalenv()))' "$config"
}
case "$mode" in
  16S|16s) render_amplicon 16S ;;
  ITS|its) render_amplicon ITS ;;
  shotgun|wgs) render_shotgun ;;
  all) render_amplicon 16S; render_amplicon ITS; render_shotgun ;;
  *) echo "Использование: ./run.sh {16S|ITS|shotgun|all} [config.yaml]" >&2; exit 2 ;;
esac
