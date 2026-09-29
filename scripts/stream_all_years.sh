#!/bin/bash
# Stream Open Payments years sequentially. API metastore exposes 2019-2025 only;
# Part D runs through 2024, so 2019-2024 is the usable overlap.
cd "$(dirname "$0")/.."
for y in "$@"; do
  if [ -s "data/raw/openpayments_diabetes_${y}.csv" ]; then
    echo "[$y] exists, skipping"; continue
  fi
  echo "=== $y ==="
  python3 scripts/stream_open_payments.py --year "$y" --class diabetes 2>&1 | grep -E 'DONE|Error|Traceback'
done
echo "ALL YEARS COMPLETE"
