#!/usr/bin/env bash
set -euo pipefail

SMR_BIN="tools/smr_linux_x86_64"
EQTL_DIR="data/eqtl/GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite"
LD_BFILE="data/reference/1000G_EUR"
GWAS_SUMMARY="data/smr/gout_ma.txt"
PROBE_LIST="data/smr/selected_probes.txt"
OUT_DIR="results/smr_gtex"

mkdir -p "${OUT_DIR}"

# Run SMR for each GTEx tissue.
for besd_file in "${EQTL_DIR}"/*.lite.besd; do
  tissue=$(basename "${besd_file}" .lite.besd)
  echo "Running SMR for ${tissue}"

  "${SMR_BIN}" \
    --bfile "${LD_BFILE}" \
    --gwas-summary "${GWAS_SUMMARY}" \
    --beqtl-summary "${EQTL_DIR}/${tissue}.lite" \
    --diff-freq-prop 0.4 \
    --peqtl-smr 5e-5 \
    --extract-probe "${PROBE_LIST}" \
    --out "${OUT_DIR}/${tissue}_gout" \
    --thread-num 32 \
    > "${OUT_DIR}/${tissue}.log" 2>&1
done

echo "SMR analysis completed."
