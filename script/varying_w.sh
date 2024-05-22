# k = 31
for w in {2..8}
do
    ./density -i random.10M.fa -k 31 -w $w -a mod-minimizer --stream 2>> mod_mini.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.fa -k 127 -w $w -a mod-minimizer --stream 2>> mod_mini.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.fa -k 31 -w $w -a lr-minimizer --stream 2>> lr_mini.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.fa -k 127 -w $w -a lr-minimizer --stream 2>> lr_mini.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.fa -k 31 -w $w -a minimizer --stream 2>> mini.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.fa -k 127 -w $w -a minimizer --stream 2>> mini.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.fa -k 31 -w $w -a miniception --stream 2>> miniception.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.fa -k 127 -w $w -a miniception --stream 2>> miniception.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.fa -k 31 -w $w -a rot-minimizer --stream 2>> rot-minimizer.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.fa -k 127 -w $w -a rot-minimizer --stream 2>> rot-minimizer.varying_w.k.127.txt
done
