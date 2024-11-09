# Sampling algorithms

This library provides a collection of sampling algorithms, including:

- random minimizers (`'M'``),
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

The code has been used for the experiments of the paper [*"The open-closed mod-minimizer algorithm"*](https://www.biorxiv.org/content/10.1101/2024.11.02.621600v1), based on the previous paper [*"The mod-minimizer: a simple and efficient sampling algorithm for long k-mers"*](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2024.11) (WABI 2024, [code](https://github.com/jermp/minimizers/releases/tag/v1.0.0)).

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
