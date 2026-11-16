#!/usr/bin/env bash
#SBATCH --job-name=vgp-defaults-smoke
#SBATCH --partition=workers
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --export=NONE
#SBATCH --output=/moosefs/erikg/vgp/logs/vgp-defaults-smoke-%j.out
#SBATCH --error=/moosefs/erikg/vgp/logs/vgp-defaults-smoke-%j.err

# PURPOSE: pure-defaults SweepGA smoke test on menidia (P08 species, 0.58 Gbp H1).
# Question (operator): are DEFAULT parameterizations reasonable, on the right
# cores, with /scratch as scratch? Runs sweepga with ONLY --num-mappings 1:1
# and --threads (aligner defaults to fastga; no --overlap/--scoring/
# --scaffold-jump/--map-pct-identity overrides). Output + intermediates live
# entirely in node-local /scratch. Comparable pair: tier3a remap production
# PAF (wfmash flags) for the same species.

set -euo pipefail

ROOT=/moosefs/erikg/lewontin-paradox
SWEEPGA=/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/bin-1/sweepga
EXPECT_SHA=fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b
H1=/moosefs/erikg/tier3data/tier3a-acquisition-20260716/menidia_menidia_fMenMen1/h1.fna
H2=/moosefs/erikg/tier3data/tier3a-acquisition-20260716/menidia_menidia_fMenMen1/h2.fna
OUT=/moosefs/erikg/tier3scratch/sweepga-defaults-smoke/menidia

[[ -x $SWEEPGA ]] || { echo "sweepga missing" >&2; exit 3; }
SHA=$(sha256sum "$SWEEPGA" | cut -d' ' -f1)
[[ $SHA == "$EXPECT_SHA" ]] || { echo "sweepga sha mismatch: $SHA" >&2; exit 3; }
[[ -s $H1 && -s $H2 ]] || { echo "inputs missing" >&2; exit 3; }
mkdir -p "$OUT"

SCRATCH=$(mktemp -d -- /scratch/vgp-defaults-smoke-${SLURM_JOB_ID}-XXXXXX)
export TMPDIR="$SCRATCH" TMP="$SCRATCH" TEMP="$SCRATCH"
trap '[[ $SCRATCH == /scratch/vgp-defaults-smoke-* ]] && rm -rf "$SCRATCH"' EXIT

cp --reflink=auto -- "$H1" "$SCRATCH/h1.fa"
cp --reflink=auto -- "$H2" "$SCRATCH/h2.fa"
cd "$SCRATCH"

echo "[$(date -u +%FT%TZ)] START defaults mapping: threads=$SLURM_CPUS_PER_TASK"
# workers nodes lack /usr/bin/time; sweepga self-reports rss/disk/peak_disk
# in its stderr, plus a light ps sampler for core utilization:
( while true; do ps -e -o pcpu=,rss=,comm= | grep -E "^( *[0-9.]+ ){2}(sweepga|fastga|wfmash)" >> "$OUT/cpu_sample.log"; sleep 5; done ) &
SAMPLER=$!
"$SWEEPGA" ./h2.fa ./h1.fa \
    --output-file "$SCRATCH/h2_to_h1.defaults.paf" \
    --num-mappings 1:1 --threads "$SLURM_CPUS_PER_TASK" \
    > "$OUT/sweepga.stdout" 2> "$OUT/sweepga.stderr"
kill $SAMPLER 2>/dev/null || true
echo "[$(date -u +%FT%TZ)] DONE mapping"

cp "$SCRATCH/h2_to_h1.defaults.paf" "$OUT/h2_to_h1.defaults.paf"
du -sb "$OUT" >> "$OUT/scratch_note.txt"
sha256sum "$OUT/h2_to_h1.defaults.paf" > "$OUT/paf.sha256"
wc -l "$OUT/h2_to_h1.defaults.paf"
echo "[$(date -u +%FT%TZ)] PROMOTED $OUT"
