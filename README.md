# Sampling algorithms

This library provides a collection of sampling algorithms, including:

- random minimizers (`'M'`),
- closed sycnmers (`'C'`, a.k.a. "miniception"),
- open syncmers (`'O'`),
- open-closed minimizers (`'OC'`),
- double decycling set based (`'DD'`),
- mod-sampling paired with any of the previous methods, thus:
	- random mod-minimizers (`'mod-M'`),
	- closed mod-minimizers (`'mod-C'`),
	- open mod-minimizers (`'mod-O'`),
	- open-closed mod-minimizers (`'mod-OC'`),
	- double-decycling mod-minimizers (`'mod-DD'`).

The library also provides the five **canonical** schemes
from the paper *"Canonical Schemes and Minimizers"*, by D. Cologni and G. E. Pibiri,
implemented in `include/canonical_schemes.hpp`:

- `'LvR'` — the minimizer of the canonical m-mers,
- `'RhL'` — the lex-smaller of the two strand minimizers (the canonical scheme of [SSHash](https://github.com/jermp/sshash)),
- `'RR'` — the minimizer under the key min{h(y), h(rc(y))},
- `'LL'` — the lex-smallest canonical m-mer,
- `'CK'` — the minimizer of the canonical k-mer.

**Note**: these schemes (and the tool `canonical_density`) use the (k, m) notation
of that paper — a window is a single k-mer, containing w = k - m + 1 m-mers — rather
than the (w, k) notation assumed elsewhere in this repository. The two conventions
match under k → m and w → k - m + 1.

The code has been used for the experiments of the paper [*"The open-closed mod-minimizer algorithm"*](https://almob.biomedcentral.com/articles/10.1186/s13015-025-00270-0) (AMB 2025), based on the previous paper [*"The mod-minimizer: a simple and efficient sampling algorithm for long k-mers"*](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2024.11) (WABI 2024, [code](https://github.com/jermp/minimizers/releases/tag/v1.0.0)).

To reproduce the experiments in the paper: first compile the code as explained below and then run the scripts [here](https://github.com/jermp/minimizers/tree/main/script#experiments).

### Compile the code

Before compiling, pull all dependencies with

	git submodule update --init --recursive

Compile the code with

    mkdir build
    cd build
    cmake ..
    make

### Quick start

After compilation, generate some random sequence (in the following example, of 1 million nucleotides) with

    ./generate_random_sequence -o test.bin -n 1000000 -s 4

and evaluate density of different methods with the tool `density`.
Some examples below.

	./density -i test.bin -k 63 -w 8 -a M --stream

	  is_forward = YES
	  num_sampled_kmers = 222230
	  num_kmers = 999938
	  num_windows = 999931
	  density = 0.222244
	  1.77795X away from lower bound 1/w = 0.125

	./density -i test.bin -k 63 -w 8 -a OC --stream

	  is_forward = YES
	  num_sampled_kmers = 181061
	  num_kmers = 999938
	  num_windows = 999931
	  density = 0.181072
	  1.44858X away from lower bound 1/w = 0.125

	./density -i test.bin -k 63 -w 8 -a mod-OC --stream

	  is_forward = YES
	  num_sampled_kmers = 138011
	  num_kmers = 999938
	  num_windows = 999931
	  density = 0.13802
	  1.10416X away from lower bound 1/w = 0.125

The density of the canonical schemes is evaluated with the tool `canonical_density`,
which takes the parameters `-k` and `-m` (see the note on notation above) and reports
the density of the scheme, whether it is forward, the number of sampled positions,
and the number of super-k-mers. For example:

	./canonical_density -i test.bin -k 31 -m 13 -a LvR

	  k = 31, m = 13, w = k - m + 1 = 19
	  is_forward = YES
	  num_sampled_positions = 100141
	  num_super_kmers = 100141
	  num_positions = 999988
	  num_windows = 999970
	  density = 0.100142 (1.00142 x d, with d = 2/(w+1) = 0.1)
	  super_kmer_rate = 0.100144 (1.00144 x d)

	./canonical_density -i test.bin -k 31 -m 13 -a RhL

	  k = 31, m = 13, w = k - m + 1 = 19
	  is_forward = NO
	  num_sampled_positions = 120840
	  num_super_kmers = 127641
	  num_positions = 999988
	  num_windows = 999970
	  density = 0.120841 (1.20841 x d, with d = 2/(w+1) = 0.1)
	  super_kmer_rate = 0.127645 (1.27645 x d)
