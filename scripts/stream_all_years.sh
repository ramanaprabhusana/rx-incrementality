#!/bin/bash
# Stream Open Payments years sequentially. API metastore exposes 2019-2025 only;
# Part D runs through 2024, so 2019-2024 is the usable overlap.
#
# A year's file exists only if its download completed (the stream writes to a
# .part file and renames it at the end), so skipping existing files is safe.
#
# Usage: scripts/stream_all_years.sh 2019 2020 ...   (CLASS=anticoagulant to stream another class)
set -o pipefail
cd "$(dirname "$0")/.." || exit 1
CLASS="${CLASS:-diabetes}"
failed=0
for y in "$@"; do
  if [ -s "data/raw/openpayments_${CLASS}_${y}.csv" ]; then
    echo "[$y] exists, skipping"; continue
  fi
  echo "=== $y ==="
  if ! python3 scripts/stream_open_payments.py --year "$y" --class "$CLASS" 2>&1 \
      | grep -E 'DONE|INCOMPLETE|NO MATCHES|Error|Traceback'; then
    echo "[$y] FAILED"; failed=1
  fi
done
if [ "$failed" -ne 0 ]; then
  echo "SOME YEARS FAILED: re-run to retry them"; exit 1
fi
echo "ALL YEARS COMPLETE"
