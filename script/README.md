# Experiments

First generate some strings of i.i.d. random characters (here, of 1 million characters from an alphabet size of 4 and 256, respectively) with

    ./generate_random_sequence -o random.1M.s4.bin -n 10000000 -s 4
    ./generate_random_sequence -o random.1M.s256.bin -n 10000000 -s 256

and then, from within the directory where the `density` tools is available (e.g., from `./build`), run

    bash ../script/varying_k.sh

## Canonical schemes

To replicate Table 3 of the paper *"Canonical Schemes and Minimizers"*, run,
from within the directory where the `canonical_density` and
`generate_random_sequence` tools are available (e.g., from `./build`),

    ../script/canonical_table3.sh

The script generates its own random strings and prints the density and the
super-k-mer rate of the five canonical schemes, averaged over independent runs.
It takes no arguments; a few environment variables tune it, e.g.

    NUM_RUNS=10 LENGTH=10000000 JOBS=8 ../script/canonical_table3.sh
    PAIRS="31:13 63:21" ../script/canonical_table3.sh
    KEEP=1 ../script/canonical_table3.sh   # also write the per-run measurements to CSV

Note that this script, like `canonical_density`, uses the (k, m) notation of
that paper: a window is a single k-mer, containing w = k - m + 1 m-mers.
