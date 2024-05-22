# Experiments

First generate a string of i.i.d. random characters (here, of 10 million characters) with

    ./generate_random_fasta -o random.10M.fa -n 10000000

and then, from within `/build`, run

    bash ../script/varying_k.sh
    bash ../script/varying_w.sh