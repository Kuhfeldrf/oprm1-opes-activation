#!/usr/bin/env bash
#SBATCH --job-name=oprm1_pack
#SBATCH --partition=normal
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/pack_%x_%j.log
# Pack one start complex into POPC:CHL1 7:3 (30 mol% cholesterol), 150 mM NaCl.
# Usage: sbatch -J pack_<start> systems/pack.sh <start> [z_box_A]
#        (active | inactive | v0.3: nb39_dam nb39_nlx free_dam free_nlx)
#
# z_box_A pins the box z axis with --dims instead of letting packmol-memgen size it.
# Needed for the Nb39 systems and for nothing else: the bilayer is pinned to
# z = [-18, +18] and the box is symmetric about it, but memgen sizes z from the
# extracellular side only (it chose +/-43 A) and ignores how far Nb39 hangs below the
# membrane. The fixed complex then lies outside its own PBC box and packmol stops with
# "Fixed molecule are outside the PBC box" (STOP 171, 13 s, no output). Nb39 reaches
# z = -54.3 A, so the half-height has to cover that plus the dist_wat layer: 150 A.
# xy stays 95 A, as --distxy_fix sets it, so lateral lipid counts stay comparable.
#
# Identical settings for both starts; --distxy_fix pins the same lateral box so the
# two systems carry (near-)identical lipid counts. Settings inherited from the parent
# repo's validated DAMGO build (docs/packing_calibration.md there):
#   --pbc mandatory (no periodic-image overlaps), apl_offset 1.15 (receptor footprint,
#   avoids over-packing), nloop_all/nloop 200/40 (NOT the 100/20 defaults).
set -euo pipefail
START=$1
ZBOX=${2:-}
if [ -n "$ZBOX" ]; then BOX=(--dims 95 95 "$ZBOX"); else BOX=(--distxy_fix 95); fi
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pilot
OUT=$ROOT/systems/$START
mkdir -p "$OUT" && cd "$OUT"
cp "$ROOT/structures/clean/$START/complex_opm.pdb" .
packmol-memgen \
  --pdb complex_opm.pdb \
  --lipids POPC:CHL1 --ratio 7:3 \
  --apl_offset 1.15 \
  "${BOX[@]}" --dist_wat 20 \
  --nloop_all 200 --nloop 40 \
  --pbc --preoriented \
  --salt --saltcon 0.15 --salt_c Na+ \
  --notprotonate --nottrim --keepligs \
  --output packed.pdb --log packmol_memgen.log
grep -q "Success" packmol.log 2>/dev/null && echo "PACKMOL: converged" \
  || echo "PACKMOL: did NOT reach tolerance (best solution written)"
for r in POPC PC PA OL CHL WAT Na+ Cl- DAM NLX; do
  printf "%-5s %s\n" "$r" "$(awk -v r=$r '$1~/ATOM|HETATM/ && $4==r' packed.pdb | awk '{print $5$6}' | sort -u | wc -l)"
done
echo "atoms $(grep -c -E '^(ATOM|HETATM)' packed.pdb)"
