#!/bin/bash
#
# Replicates Table 3 of the paper
#
#     D. Cologni and G. E. Pibiri, "Canonical Schemes and Minimizers".
#
# For each (k, m) pair it runs the five canonical schemes -- LvR, RhL, RR, LL
# and CK -- over NUM_RUNS independent random DNA strings of LENGTH symbols,
# and reports the mean of each measured quantity, together with the standard
# deviation over the runs.
#
# NOTE ON NOTATION: the (k, m) notation of the paper is used, hence -k and -m;
# a window is a single k-mer, containing w = k - m + 1 m-mers. The rest of this
# repository uses the (w, k) notation, matching under k -> m, w -> k - m + 1.
#
# Run it from the build directory, e.g.
#
#     cd build && ../script/canonical_table3.sh
#
# Environment variables (all optional):
#
#     NUM_RUNS   number of independent random strings          (default 5)
#     LENGTH     length of each random string                  (default 1000000)
#     JOBS       number of measurements to run concurrently    (default: nproc)
#     SEED0      first sequence seed; run r uses SEED0 + r     (default 1)
#     PAIRS      space-separated "k:m" list, overriding the default set
#     KEEP       if non-empty, keep the generated sequences and the raw CSV
#
set -eu

NUM_RUNS=${NUM_RUNS:-5}
LENGTH=${LENGTH:-1000000}
JOBS=${JOBS:-$( (command -v nproc >/dev/null && nproc) || echo 4)}
SEED0=${SEED0:-1}

# The (k, m) pairs. Both odd and even k are included: the tie-break of
# Proposition 24 is mirror-equivariant only for odd k, and the table shows that
# the choice of parity leaves the densities unaffected. Every m clears the
# genericity condition (3) of the paper, m >= (4 + eps) log_sigma w, which at
# sigma = 4 reads m >= (2 + eps/2) log_2 w.
PAIRS=${PAIRS:-"16:9 17:9 30:13 31:13 63:21 64:21 100:25 102:25 127:25 128:25"}

SCHEMES="LvR RhL RR LL CK"

for exe in ./canonical_density ./generate_random_sequence; do
    if [ ! -x "$exe" ]; then
        echo "error: $exe not found; run this script from the build directory" >&2
        exit 1
    fi
done

workdir=$(mktemp -d "${TMPDIR:-/tmp}/canonical_table3.XXXXXX")
csv="$workdir/measurements.csv"
cleanup() { if [ -z "${KEEP:-}" ]; then rm -rf "$workdir"; fi; }
trap cleanup EXIT

echo "# Table 3: $NUM_RUNS random strings of $LENGTH symbols, sigma = 4, $JOBS jobs" >&2

# 1. Generate the independent random strings, one per run.
for r in $(seq 1 "$NUM_RUNS"); do
    ./generate_random_sequence -o "$workdir/seq.$r.bin" -n "$LENGTH" -s 4 \
                               --seed $((SEED0 + r)) >/dev/null
done

# 2. Measure. One job per (pair, scheme, run); each writes a single CSV line,
#    to its own file so that concurrent jobs cannot interleave their output.
measure() {
    k=${1%%:*}; m=${1##*:}; alg=$2; r=$3; workdir=$4
    ./canonical_density -i "$workdir/seq.$r.bin" -k "$k" -m "$m" -a "$alg" | awk \
        -v k="$k" -v m="$m" -v alg="$alg" '
        /^  is_forward =/      { fwd = ($3 == "YES") ? 1 : 0 }
        /^  num_super_kmers =/ { skm_abs = $3 }
        /^  density =/         { den = $3; dr  = $4; gsub(/[()]/, "", dr) }
        /^  super_kmer_rate =/ { skm = $3; skr = $4; gsub(/[()]/, "", skr) }
        END { printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
                     k, m, alg, den, dr, skm, skr, skm_abs, fwd }
    ' > "$workdir/out.$k.$m.$alg.$r.csv"
}
export -f measure

for pair in $PAIRS; do
    for alg in $SCHEMES; do
        for r in $(seq 1 "$NUM_RUNS"); do echo "$pair $alg $r $workdir"; done
    done
done | xargs -P "$JOBS" -L 1 bash -c 'measure $0 $1 $2 $3'

cat "$workdir"/out.*.csv > "$csv"

# 3. Aggregate: mean and standard deviation over the runs, then lay the table
#    out in the columns of Table 3. LvR, RR and LL are forward, so for them
#    density and super-k-mer rate are equal and only the density is shown; the
#    super-k-mer count of CK is given in absolute terms, being Theta(1) rather
#    than a multiple of d.
awk -F, -v runs="$NUM_RUNS" '
    function mean(s, n)  { return s / n }
    function sdev(s, ss, n,   v) { v = ss / n - (s / n) ^ 2; return v > 0 ? sqrt(v) : 0 }
    function acc(key, x) { sum[key] += x; sqsum[key] += x * x; cnt[key] += 1 }
    {
        k = $1; m = $2; alg = $3
        key = k "," m
        if (!(key in seen)) { seen[key] = 1; order[++nkeys] = key }
        acc(key "," alg ",den",  $4)   # density
        acc(key "," alg ",dr",   $5)   # density / d
        acc(key "," alg ",skm",  $6)   # super-k-mer rate
        acc(key "," alg ",skr",  $7)   # super-k-mer rate / d
        acc(key "," alg ",abs",  $8)   # super-k-mer rate, absolute
        acc(key "," alg ",fwd",  $9)   # forwardness
        if ($9 == 0) notfwd[key "," alg] = 1
    }
    END {
        # Order the rows by increasing w, as in the paper, rather than by the
        # lexicographic order in which the per-job files happened to be read.
        for (i = 2; i <= nkeys; ++i) {
            key = order[i]; split(key, kv, ","); wi = kv[1] - kv[2] + 1
            for (j = i - 1; j >= 1; --j) {
                split(order[j], kv, ","); if (kv[1] - kv[2] + 1 <= wi) break
                order[j + 1] = order[j]
            }
            order[j + 1] = key
        }
        printf "\n"
        printf "%-4s %-3s %-4s %-8s | %-15s | %-15s | %-6s | %-17s\n",
               "", "", "", "", "RhL", "CK", "LL", "LvR / RR"
        printf "%-4s %-3s %-4s %-8s | %-7s %-7s | %-7s %-7s | %-6s | %-8s %-8s\n",
               "k", "m", "w", "d", "d/d", "skm/d", "d/d", "skm", "d/d", "d(LvR)", "d(RR)"
        printf "%s\n", "---------------------------------------------------------------------------------"
        for (i = 1; i <= nkeys; ++i) {
            key = order[i]
            split(key, kv, ",")
            k = kv[1]; m = kv[2]; w = k - m + 1
            d = 2.0 / (w + 1)
            printf "%-4s %-3s %-4d %-8.5f | %-7.3f %-7.3f | %-7.3f %-7.5f | %-6.3f | %-8.5f %-8.5f\n",
                   k, m, w, d,
                   mean(sum[key ",RhL,dr"], runs), mean(sum[key ",RhL,skr"], runs),
                   mean(sum[key ",CK,dr"],  runs), mean(sum[key ",CK,skm"],  runs),
                   mean(sum[key ",LL,dr"],  runs),
                   mean(sum[key ",LvR,den"], runs), mean(sum[key ",RR,den"], runs)
        }
        printf "%-4s %-3s %-4s %-8s | %-7s %-7s | %-7s %-7s | %-6s | %-8s %-8s\n",
               "thry", "", "", "", "", "-> 4/3", "-> 2", "", "", "d", "d"

        # The ratio the theory pins down: skm(RhL) / skm(LvR) -> 4/3.
        printf "\n%-4s %-3s %-4s  %-19s  %s\n",
               "k", "m", "w", "skm(RhL)/skm(LvR)", "max sdev of d/d over the " runs " runs"
        printf "%s\n", "---------------------------------------------------------------------"
        for (i = 1; i <= nkeys; ++i) {
            key = order[i]
            split(key, kv, ",")
            k = kv[1]; m = kv[2]; w = k - m + 1
            ratio = mean(sum[key ",RhL,skm"], runs) / mean(sum[key ",LvR,skm"], runs)
            worst = 0
            for (a = 1; a <= 5; ++a) {
                alg = (a == 1 ? "LvR" : a == 2 ? "RhL" : a == 3 ? "RR" : a == 4 ? "LL" : "CK")
                s = sdev(sum[key "," alg ",dr"], sqsum[key "," alg ",dr"], runs)
                if (s > worst) worst = s
            }
            printf "%-4s %-3s %-4d  %-19.4f  %.5f\n", k, m, w, ratio, worst
        }

        # Forwardness: LvR, RR and LL must be forward; RhL and CK must not be.
        printf "\nforwardness (expected: LvR, RR, LL forward; RhL, CK not)\n"
        bad = 0
        for (i = 1; i <= nkeys; ++i) {
            key = order[i]
            for (a = 1; a <= 5; ++a) {
                alg = (a == 1 ? "LvR" : a == 2 ? "RhL" : a == 3 ? "RR" : a == 4 ? "LL" : "CK")
                expect = (alg == "RhL" || alg == "CK") ? 0 : 1
                got = mean(sum[key "," alg ",fwd"], runs)
                if (got != expect) {
                    printf "  UNEXPECTED: k,m = %s scheme %s: forward on %g of %d runs\n",
                           key, alg, got * runs, runs
                    bad = 1
                }
            }
        }
        if (!bad) printf "  as expected on every (k, m) and every run\n"
    }
' "$csv"

if [ -n "${KEEP:-}" ]; then
    cp "$csv" ./canonical_table3.csv
    echo "" >&2
    echo "# raw per-run measurements written to ./canonical_table3.csv" >&2
fi
