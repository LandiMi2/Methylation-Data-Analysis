#using deeptools - https://deeptools.readthedocs.io/en/latest/index.html

#global methylation profiles 
multiBigwigSummary bins -b *.bw --binSize 100000 -out meth100kbCpG.npz --outRawCounts meth100kbCpG.tab
#meth100kbCpG.tab loaded in R to plot krayoploteR

#generate metaplots for genes 
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

##PCA plots 
multiBigwigSummary bins \
  -b CMV1.CpG.bw CMV2.CpG.bw CaMV1.CpG.bw CaMV2.CpG.bw Mock1.CpG.bw Mock2.CpG.bw \
  --binSize 10000 \
  -o CpG.globalMethylation.npz \
  --outRawCounts CpG.globalMethylation.tab

plotPCA \
  -in CpG.globalMethylation.npz \
  -o CpG.methylation.global.PCA.pdf \
  --plotTitle "PCA" \
  --labels "CMV1" "CMV2" "CaMV1" "CaMV2" "Mock1" "Mock2" \
  --colors "#ff2a0e" "#ff2a0e" "#ff9f4e" "#ff9f4e" "#1f77b4" "#1f77b4"

##TE metaplot
computeMatrix scale-regions -S ../CMV1.CpG.bw ../CMV2.CpG.bw ../CMV3.CpG.bw \
../CaMV1.CpG.bw ../CaMV2.CpG.bw ../CaMV3.CpG.bw \
../Mock1.CpG.bw ../Mock2.CpG.bw ../Mock3.CpG.bw \
-R te.bed --beforeRegionStartLength 2000 -p 10 --afterRegionStartLength 2000 \
--regionBodyLength 5000 --samplesLabel "CMV1" "CMV2" "CMV3" "CaMV1" "CaMV2" "CaMV3" "Mock1" "Mock2" "Mock3" -o CpG.te.matrix.gz

#metaplot
plotProfile -m  CpG.te.matrix.gz -out CpG.TE.pdf --perGroup \
--colors "#ff2a0e" "#ff2a0e" "#ff2a0e" "#ff9f4e" "#ff9f4e" "#ff9f4e" "#1f77b4" "#1f77b4" "#1f77b4" \
--plotTitle "Transposon" \
--refPointLabel "TSS/TTS"

