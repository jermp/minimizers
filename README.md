# Sampling algorithms

This library provides a collection of sampling algorithms, including:

- mod-sampling,
- lr-minimizers (a "context-sensitive" version of closed syncmers),
- mod-minimizers,
- "classic" minimizers,
- miniception,
- rotational-minimizers,
- decycling set based minimizers,
- closed-syncmers,
- open-syncmers,
- open-closed-syncmers.

The code has been used for the experiments of the paper [*"The mod-minimizer: a simple and efficient sampling algorithm for long k-mers"*](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2024.11), published in WABI 2024.

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

	./density -i test.bin -k 63 -w 8 -a minimizer --stream

	  num_sampled_kmers = 222287
	  num_kmers = 999938
	  num_windows = 999931
	  density = 0.222301
	  1.77841X away from lower bound 1/w = 0.125
	  calculation using closed formulas:
	     density = 0.222222
	     1.77778X away from lower bound 1/w = 0.125

	./density -i test.bin -k 63 -w 8 -a lr-minimizer --stream

	  num_sampled_kmers = 176558
	  num_kmers = 999938
	  num_windows = 999931
	  density = 0.176569
	  1.41255X away from lower bound 1/w = 0.125
	  calculation using closed formulas:
	     density = 0.176471
	     1.41176X away from lower bound 1/w = 0.125

	./density -i test.bin -k 63 -w 8 -a mod-minimizer --stream

	  num_sampled_kmers = 138501
	  num_kmers = 999938
	  num_windows = 999931
	  density = 0.13851
	  1.10808X away from lower bound 1/w = 0.125
	  calculation using closed formulas:
	     density = 0.138462
	     1.10769X away from lower bound 1/w = 0.125