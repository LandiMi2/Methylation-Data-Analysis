#filter 5 cov
for file in *.cov; do
    base=$(basename "$file" .cov)
    awk '($5 + $6) >= 5 {
        $2 = $2 - 1;
        $4 = sprintf("%.0f", $4);
        print
    }' OFS="\t" "$file" > "${base}.bedGraph"
done

#sort
for file in *.bedGraph; do
    base=$(basename "$file" .bedGraph)

    #sort and extract
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' "$file" | sort -k1,1 -k2,2n > "${base}.CpG.bedGraph"
done

#common sites
bedtools unionbedg -i Mock1.CpG.bedGraph Mock2.CpG.bedGraph Mock3.CpG.bedGraph \
-filler - -names Mock1 Mock2 Mock3 > unionCpG.tsv


#remove
grep -v $'\t-' unionCpG.tsv > unionCpG_common.tsv 

cut -f1-3,4 unionCpG_common.tsv > Mock1_common_CpG.txt
cut -f1-3,5 unionCpG_common.tsv > Mock2_common_CpG.txt
cut -f1-3,6 unionCpG_common.tsv > Mock3_common_CpG.txt

#get matrix file 
for file in *_common_CpG.txt; do
    base=$(basename "$file" _common_CpG.txt)
    #convert to bigwig
    bedGraphToBigWig "$file" chrom.sizes  "${base}.CpG.bw"
done
