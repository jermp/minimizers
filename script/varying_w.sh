# k = 31
for w in {2..8}
do
    ./density -i random.10M.s4.bin -k 31 -w $w -a mod-minimizer --stream 2>> mod_mini.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.s4.bin -k 127 -w $w -a mod-minimizer --stream 2>> mod_mini.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.s4.bin -k 31 -w $w -a lr-minimizer --stream 2>> lr_mini.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.s4.bin -k 127 -w $w -a lr-minimizer --stream 2>> lr_mini.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.s4.bin -k 31 -w $w -a minimizer --stream 2>> mini.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.s4.bin -k 127 -w $w -a minimizer --stream 2>> mini.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.s4.bin -k 31 -w $w -a miniception --stream 2>> miniception.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.s4.bin -k 127 -w $w -a miniception --stream 2>> miniception.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.s4.bin -k 31 -w $w -a rot-minimizer-alt --stream 2>> rot-minimizer-alt.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.s4.bin -k 127 -w $w -a rot-minimizer-alt --stream 2>> rot-minimizer-alt.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.s4.bin -k 31 -w $w -a rot-minimizer-orig --stream 2>> rot-minimizer-orig.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.s4.bin -k 127 -w $w -a rot-minimizer-orig --stream 2>> rot-minimizer-orig.varying_w.k.127.txt
done

# k = 31
for w in {2..8}
do
    ./density -i random.10M.s4.bin -k 31 -w $w -a decycling --stream 2>> decycling.varying_w.k.31.txt
done

# k = 127
for w in {2..24}
do
    ./density -i random.10M.s4.bin -k 127 -w $w -a double-decycling --stream 2>> double-decycling.varying_w.k.127.txt
done