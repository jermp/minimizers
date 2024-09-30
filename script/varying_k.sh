# w = 8
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 8 -a M --stream 2>> M.varying_k.w.8.txt
done

# w = 8
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 8 -a C --stream 2>> C.varying_k.w.8.txt
done

# w = 8
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 8 -a O --stream 2>> O.varying_k.w.8.txt
done

# w = 8
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 8 -a OC --stream 2>> OC.varying_k.w.8.txt
done

# w = 8
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 8 -a mod-M --stream 2>> mod-M.varying_k.w.8.txt
done

# w = 8
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 8 -a mod-C --stream 2>> mod-C.varying_k.w.8.txt
done

# w = 8
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 8 -a mod-O --stream 2>> mod-O.varying_k.w.8.txt
done

# w = 8
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 8 -a mod-OC --stream 2>> mod-OC.varying_k.w.8.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a M --stream 2>> M.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a C --stream 2>> C.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a O --stream 2>> O.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a OC --stream 2>> OC.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-M --stream 2>> mod-M.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-C --stream 2>> mod-C.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-O --stream 2>> mod-O.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-OC --stream 2>> mod-OC.varying_k.w.24.txt
done