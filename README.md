# Sampling algorithms

Before compiling, pull all dependencies with

	git submodule update --init --recursive

Compile the code with

    mkdir build
    cd build
    cmake ..
    make

and generate some random sequence (in the following example, of 1 million nucleotides) with

    ./generate_random_fasta -o test.fa -n 1000000

Lastly, evaluate density of different methods with the tool `density`.
Some examples below.

	./density -i test.fa -k 63 -w 8 -a minimizer --stream

	  is_forward = YES (expected = YES)
	  num_sampled_kmers = 222355
	  num_kmers = 999938
	  density = 0.222369
	  1.77895X away from lower bound 1/w = 0.125
	  calculation using closed formulas:
	     density = 0.222222
	     1.77778X away from lower bound 1/w = 0.125

	./density -i test.fa -k 63 -w 8 -a lr-minimizer --stream

	  is_forward = YES (expected = YES)
	  num_sampled_kmers = 176441
	  num_kmers = 999938
	  density = 0.176452
	  1.41162X away from lower bound 1/w = 0.125
	  calculation using closed formulas:
	     density = 0.176471
	     1.41176X away from lower bound 1/w = 0.125

	./density -i test.fa -k 63 -w 8 -a mod-minimizer --stream

	  is_forward = YES (expected = YES)
	  num_sampled_kmers = 138474
	  num_kmers = 999938
	  density = 0.138483
	  1.10786X away from lower bound 1/w = 0.125
	  calculation using closed formulas:
	     density = 0.138462
	     1.10769X away from lower bound 1/w = 0.125
