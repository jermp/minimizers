# k = 31
for w in {2..31}
do
    ./density -i test.fa -k 31 -w $w -a mod-minimizer 2>> mod_mini.varying_w.k.31.txt
done

# k = 63
for w in {2..63}
do
    ./density -i test.fa -k 63 -w $w -a mod-minimizer 2>> mod_mini.varying_w.k.63.txt
done

# k = 31
for w in {2..31}
do
    ./density -i test.fa -k 31 -w $w -a lr-minimizer 2>> lr_mini.varying_w.k.31.txt
done

# k = 63
for w in {2..63}
do
    ./density -i test.fa -k 63 -w $w -a lr-minimizer 2>> lr_mini.varying_w.k.63.txt
done

# k = 31
for w in {2..31}
do
    ./density -i test.fa -k 31 -w $w -a minimizer 2>> mini.varying_w.k.31.txt
done

# k = 63
for w in {2..63}
do
    ./density -i test.fa -k 63 -w $w -a minimizer 2>> mini.varying_w.k.63.txt
done

# k = 31
for w in {2..31}
do
    ./density -i test.fa -k 31 -w $w -a miniception 2>> miniception.varying_w.k.31.txt
done

# k = 63
for w in {2..63}
do
    ./density -i test.fa -k 63 -w $w -a miniception 2>> miniception.varying_w.k.63.txt
done

# k = 31
for w in {2..31}
do
    ./density -i test.fa -k 31 -w $w -a rot-minimizer 2>> rot-minimizer.varying_w.k.31.txt
done

# k = 63
for w in {2..63}
do
    ./density -i test.fa -k 63 -w $w -a rot-minimizer 2>> rot-minimizer.varying_w.k.63.txt
done