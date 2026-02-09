#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${HERE}/data"
mkdir -p "${DATA_DIR}"

# Pin to a specific upstream commit for reproducibility.
BESTEST_GSR_COMMIT="c82e3718e256e18e404629807b4398b42b7e1639"

EPW_NAME="USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw"
EPW_URL="https://raw.githubusercontent.com/NREL/BESTEST-GSR/${BESTEST_GSR_COMMIT}/integration_testing/weather/${EPW_NAME}"

if [[ -f "${DATA_DIR}/${EPW_NAME}" ]]; then
  echo "Already present: ${DATA_DIR}/${EPW_NAME}"
  exit 0
fi

echo "Downloading EPW:"
echo "  ${EPW_URL}"

curl -L --fail -o "${DATA_DIR}/${EPW_NAME}" "${EPW_URL}"

echo "Saved:"
echo "  ${DATA_DIR}/${EPW_NAME}"

