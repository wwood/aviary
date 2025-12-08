#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "CONDA_PREFIX is not set; unable to determine pixi environment location" >&2
  exit 1
fi

GATB_DIR="${CONDA_PREFIX}/share/gatb-minia-pipeline"
GATB_BIN="${CONDA_PREFIX}/bin/gatb"

if [[ -x "${GATB_BIN}" ]]; then
  exit 0
fi

rm -rf "${GATB_DIR}"
mkdir -p "${GATB_DIR}"

GIT_SSL_NO_VERIFY=1 git clone --depth 1 https://github.com/GATB/gatb-minia-pipeline.git "${GATB_DIR}" >/dev/null
ln -sf "${GATB_DIR}/gatb" "${GATB_BIN}"
chmod +x "${GATB_BIN}"
