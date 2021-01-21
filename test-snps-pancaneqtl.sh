# run this script on the Git bash terminal and not on cygwin!
# there are no SNPs for hair color overlapping trans-eQTLs.

###############################################################
# eqtl="cis"      # options: "trans" "cis"
# type="UVM" 	# options: "UVM" "SKCM"
# color="blonde"  # options: "brown" "red" "blonde"
###############################################################

### for hair color SNPs:

eqtl="cis"

for color in blonde red brown
do
	for type in UVM SKCM
	do
		cp label.txt hair-${color}-${eqtl}-eqtl-${type}.out
		lines=$(wc -l top-hair-${color}.txt | awk '{print $1}')
		for (( c=1; c<=$lines; c++ ))
		do
			snp="$(sed -n "${c}p" top-hair-${color}.txt)"
			grep -w "${snp}" ${eqtl}_eQTLs_${type}.txt >> hair-${color}-${eqtl}-eqtl-${type}.out
		done
	done
done


################################################################
###############################################################
eqtl="cis"      # options: "trans" "cis"
type="UVM" 	# options: "UVM" "SKCM"
color="blonde"  # options: "green" "hazel" "blue"
###############################################################

### for eye color SNPs:

for color in blue green hazel
do
	for eqtl in cis trans
	do
		for type in UVM SKCM
		do
			cp label.txt eye-${color}-${eqtl}-eqtl-${type}.out
			lines=$(wc -l top-eye-${color}.txt | awk '{print $1}')
			for (( c=1; c<=$lines; c++ ))
			do
				snp="$(sed -n "${c}p" top-eye-${color}.txt)"
				grep -w "${snp}" ${eqtl}_eQTLs_${type}.txt >> eye-${color}-${eqtl}-eqtl-${type}.out
			done
		done
	done
done

