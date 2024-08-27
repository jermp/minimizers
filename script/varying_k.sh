# w = 8
for k in {5..100}
do
	./density -i random.10M.fa -k $k -w 8 -a mod-minimizer --stream 2>> mod_mini.varying_k.w.8.txt
done

# w = 8
for k in {5..100}
do
	./density -i random.10M.fa -k $k -w 8 -a lr-minimizer --stream 2>> lr_mini.varying_k.w.8.txt
done

# w = 8
for k in {5..100}
do
	./density -i random.10M.fa -k $k -w 8 -a minimizer --stream 2>> mini.varying_k.w.8.txt
done

# w = 8
for k in {5..100}
do
	./density -i random.10M.fa -k $k -w 8 -a miniception --stream 2>> miniception.varying_k.w.8.txt
done

# w = 8
for k in {5..100}
do
	./density -i random.10M.fa -k $k -w 8 -a open-closed-syncmer --stream 2>> open-closed-syncmer.varying_k.w.8.txt
done

# w = 8
for k in {5..100}
do
    ./density -i random.10M.fa -k $k -w 8 -a rot-minimizer-alt --stream 2>> rot-minimizer-alt.varying_k.w.8.txt
done

# w = 8
for k in {5..100}
do
    ./density -i random.10M.fa -k $k -w 8 -a rot-minimizer-orig --stream 2>> rot-minimizer-orig.varying_k.w.8.txt
done

# w = 8
for k in {5..100}
do
	./density -i random.10M.fa -k $k -w 8 -a decycling --stream 2>> decycling.varying_k.w.8.txt
done

# w = 8
for k in {5..100}
do
	./density -i random.10M.fa -k $k -w 8 -a double-decycling --stream 2>> double-decycling.varying_k.w.8.txt
done

# w = 24
for k in {5..300}
do
	./density -i random.10M.fa -k $k -w 24 -a mod-minimizer --stream 2>> mod_mini.varying_k.w.24.txt
done

# w = 24
for k in {5..300}
do
	./density -i random.10M.fa -k $k -w 24 -a lr-minimizer --stream 2>> lr_mini.varying_k.w.24.txt
done

# w = 24
for k in {5..300}
do
	./density -i random.10M.fa -k $k -w 24 -a minimizer --stream 2>> mini.varying_k.w.24.txt
done

# w = 24
for k in {5..300}
do
	./density -i random.10M.fa -k $k -w 24 -a miniception --stream 2>> miniception.varying_k.w.24.txt
done

# w = 24
for k in {5..300}
do
	./density -i random.10M.fa -k $k -w 24 -a open-closed-syncmer --stream 2>> open-closed-syncmer.varying_k.w.24.txt
done

# w = 24
for k in {5..300}
do
    ./density -i random.10M.fa -k $k -w 24 -a rot-minimizer-alt --stream 2>> rot-minimizer-alt.varying_k.w.24.txt
done

# w = 24
for k in {5..300}
do
    ./density -i random.10M.fa -k $k -w 24 -a rot-minimizer-orig --stream 2>> rot-minimizer-orig.varying_k.w.24.txt
done

# w = 24
for k in {5..300}
do
	./density -i random.10M.fa -k $k -w 24 -a decycling --stream 2>> decycling.varying_k.w.24.txt
done

# w = 24
for k in {5..300}
do
	./density -i random.10M.fa -k $k -w 24 -a double-decycling --stream 2>> double-decycling.varying_k.w.24.txt
done
