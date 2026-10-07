#!/usr/bin/env bash
# Download the runbook §4 entries (mmCIF) from RCSB into structures/raw/.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p raw
for id in 9MQH 9MQJ 9MQI 9PXU 10TM 9PPQ 8F7Q 8Y72; do
  [ -s "raw/$id.cif" ] || curl -sSfL -o "raw/$id.cif" "https://files.rcsb.org/download/$id.cif"
  echo "$id $(wc -c < raw/$id.cif) bytes"
done
