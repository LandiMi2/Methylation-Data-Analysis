computeMatrix scale-regions -S CMV1.CpG.bw CMV2.CpG.bw CMV3.CpG.bw \
CaMV1.CpG.bw CaMV2.CpG.bw CaMV3.CpG.bw \
Mock1.CpG.bw Mock2.CpG.bw Mock3.CpG.bw \
-R genes.bed --beforeRegionStartLength 2000 -p 10 --afterRegionStartLength 2000 \
--regionBodyLength 5000 --samplesLabel "CMV1" "CMV2" "CMV3" "CaMV1" "CaMV2" "CaMV3" "Mock1" "Mock2" "Mock3" -o CpG.matrix.gz

#metaplot
plotProfile -m  CpG.matrix.gz -out CpG.Genes.pdf --perGroup \
--colors "#ff2a0e" "#b21807" "#67000d" "#ffb347" "#ff7f0e" "#a85200" "#1f77b4" "#084594" "#08306b" \
--plotTitle "Gene Body" \
--refPointLabel "TSS/TTS" 
