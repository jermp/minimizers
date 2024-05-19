# w = 8
for k in {5..100}
do
	./density -i test.fa -k $k -w 8 -a mod-minimizer 2>> mod_mini.varying_k.w.8.txt
done

# w = 16
for k in {5..100}
do
	./density -i test.fa -k $k -w 16 -a mod-minimizer 2>> mod_mini.varying_k.w.16.txt
done

# w = 24
for k in {5..100}
do
	./density -i test.fa -k $k -w 24 -a mod-minimizer 2>> mod_mini.varying_k.w.24.txt
done

# w = 8
for k in {5..100}
do
	./density -i test.fa -k $k -w 8 -a lr-minimizer 2>> lr_mini.varying_k.w.8.txt
done

# w = 16
for k in {5..100}
do
	./density -i test.fa -k $k -w 16 -a lr-minimizer 2>> lr_mini.varying_k.w.16.txt
done

# w = 24
for k in {5..100}
do
	./density -i test.fa -k $k -w 24 -a lr-minimizer 2>> lr_mini.varying_k.w.24.txt
done

# w = 8
for k in {5..100}
do
	./density -i test.fa -k $k -w 8 -a minimizer 2>> mini.varying_k.w.8.txt
done

# w = 16
for k in {5..100}
do
	./density -i test.fa -k $k -w 16 -a minimizer 2>> mini.varying_k.w.16.txt
done

# w = 24
for k in {5..100}
do
	./density -i test.fa -k $k -w 24 -a minimizer 2>> mini.varying_k.w.24.txt
done

# w = 8
for k in {5..100}
do
	./density -i test.fa -k $k -w 8 -a miniception 2>> miniception.varying_k.w.8.txt
done

# w = 16
for k in {5..100}
do
	./density -i test.fa -k $k -w 16 -a miniception 2>> miniception.varying_k.w.16.txt
done

# w = 24
for k in {5..100}
do
	./density -i test.fa -k $k -w 24 -a miniception 2>> miniception.varying_k.w.24.txt
done

# w = 8
for k in {5..100}
do
    ./density -i test.fa -k $k -w 8 -a rot-minimizer 2>> rot-minimizer.varying_k.w.8.txt
done

# w = 16
for k in {5..100}
do
    ./density -i test.fa -k $k -w 16 -a rot-minimizer 2>> rot-minimizer.varying_k.w.16.txt
done

# w = 24
for k in {5..100}
do
    ./density -i test.fa -k $k -w 24 -a rot-minimizer 2>> rot-minimizer.varying_k.w.24.txt
done