# Experiments

First generate a string of i.i.d. random characters (here, of 10 million characters from al alphabet size of 4, e.g. DNA) with

    ./generate_random_sequence -o random.10M.s4.bin -n 10000000 -s 4

and then, from within the directory where the `density` tools is available (e.g., `/build`), run

    bash ../script/varying_k.sh
    bash ../script/varying_w.sh