#!/bin/bash
set -euo pipefail

# ============================
# 用法说明
# ============================
usage() {
cat << EOF

Usage:
  ./install_r_packages.sh [options]

Options:
  --cran     "pkg1,pkg2,pkg3"
  --bioc     "pkg1,pkg2"
  --github   "user/repo,user2/repo2"
  --mirror   CRAN mirror (default: TUNA)

Example:
  ./install_r_packages.sh \
    --cran "optparse,dplyr,ggplot2" \
    --bioc "phyloseq,ANCOMBC" \
    --github "microbiota/amplicon"

EOF
exit 1
}

# ============================
# 默认参数
# ============================
CRAN_MIRROR="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
CRAN_PKGS=""
BIOC_PKGS=""
GITHUB_PKGS=""

# ============================
# 参数解析
# ============================
while [[ $# -gt 0 ]]; do
  case $1 in
    --cran)
      CRAN_PKGS="$2"
      shift 2
      ;;
    --bioc)
      BIOC_PKGS="$2"
      shift 2
      ;;
    --github)
      GITHUB_PKGS="$2"
      shift 2
      ;;
    --mirror)
      CRAN_MIRROR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

# ============================
# 执行 R 安装逻辑
# ============================
Rscript - << EOF
options(repos = c(CRAN="$CRAN_MIRROR"))

installed <- rownames(installed.packages())

# ---------- CRAN ----------
if ("$CRAN_PKGS" != "") {
  cran_pkgs <- unlist(strsplit("$CRAN_PKGS", ","))
  for (pkg in cran_pkgs) {
    if (!pkg %in% installed) {
      message("Installing CRAN package: ", pkg)
      install.packages(pkg)
    }
  }
}

# ---------- Bioconductor ----------
if ("$BIOC_PKGS" != "") {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  bioc_pkgs <- unlist(strsplit("$BIOC_PKGS", ","))
  for (pkg in bioc_pkgs) {
    if (!pkg %in% installed) {
      message("Installing Bioconductor package: ", pkg)
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }
}

# ---------- GitHub ----------
if ("$GITHUB_PKGS" != "") {
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

  github_pkgs <- unlist(strsplit("$GITHUB_PKGS", ","))
  for (repo in github_pkgs) {
    pkg <- sub(".*/", "", repo)
    if (!pkg %in% installed) {
      message("Installing GitHub package: ", repo)
      devtools::install_github(repo, upgrade = "never")
    }
  }
}
EOF

echo "✅ R packages installation finished"
