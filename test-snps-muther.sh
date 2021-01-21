### for hair color SNPs:

for color in blonde red brown
do
	cp label-muther.txt hair-${color}-MuTHER.out
	lines=$(wc -l top-hair-${color}.txt | awk '{print $1}')
	for (( c=1; c<=$lines; c++ ))
	do
		snp="$(sed -n "${c}p" top-hair-${color}.txt)"
		grep -w "${snp}" MuTHER_top_cis_eQTL_per_probe_Skin.txt >> hair-${color}-MuTHER.out
	done
done

### for eye color SNPs:

for color in blue green hazel
do
	cp label-muther.txt eye-${color}-MuTHER.out
	lines=$(wc -l top-eye-${color}.txt | awk '{print $1}')
	for (( c=1; c<=$lines; c++ ))
	do
		snp="$(sed -n "${c}p" top-eye-${color}.txt)"
		grep -w "${snp}" MuTHER_top_cis_eQTL_per_probe_Skin.txt >> eye-${color}-MuTHER.out
	done
done