# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a M --stream 2>> M.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a C --stream 2>> C.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a O --stream 2>> O.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a OC --stream 2>> OC.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a DD --stream 2>> DD.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-M --stream 2>> mod-M.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-C --stream 2>> mod-C.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-O --stream 2>> mod-O.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-OC --stream 2>> mod-OC.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s4.bin -k $k -w 24 -a mod-DD --stream 2>> mod-DD.s4.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a M --stream 2>> M.s256.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a C --stream 2>> C.s256.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a O --stream 2>> O.s256.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a OC --stream 2>> OC.s256.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a DD --stream 2>> DD.s256.varying_k.w.24.txt
done

# # w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a mod-M --stream 2>> mod-M.s256.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a mod-C --stream 2>> mod-C.s256.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a mod-O --stream 2>> mod-O.s256.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a mod-OC --stream 2>> mod-OC.s256.varying_k.w.24.txt
done

# w = 24
for k in {5..75}
do
	./density -i random.1M.s256.bin -k $k -w 24 -a mod-DD --stream 2>> mod-DD.s256.varying_k.w.24.txt
done
