# Experiments

First generate some strings of i.i.d. random characters (here, of 1 million characters from an alphabet size of 4 and 256, respectively) with

    ./generate_random_sequence -o random.1M.s4.bin -n 10000000 -s 4
    ./generate_random_sequence -o random.1M.s256.bin -n 10000000 -s 256

and then, from within the directory where the `density` tools is available (e.g., from `./build`), run

    bash ../script/varying_k.sh
