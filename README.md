# DNA methylation analysis 

This repository contains workflow and scripts for processing whole-genome bisulfite sequencing (WGBS) data, performing global methylation profiling, and identifying differentially methylated regions (DMRs) in Arabidopsis thaliana using the TAIR10.1 genome.  <br>

Raw bisulfite sequencing reads were assessed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and adapter trimming and quality filtering were performed using [Trim Galore](https://github.com/FelixKrueger/TrimGalore).  

Align reads to the TAIR10.1 reference genome using Bismark: <br>
Mapping using [Snakefile](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/MethPipeline/Snakefile) <br>
Generate bedGraph and coverage files using [bedGraphSnakefile](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/MethPipeline/bedGraphSnakefile)

### Global methylation profiles
Prefiltering step, getting common sites, and converting to bigWig files - [filter.sh](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/GlobalProfile/filter.sh) <br>
Generate metaplots over genes/TEs -  [Plots.sh](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/GlobalProfile/Plots.sh) <br>
Global methylation profiles using _karyoploteR_ - [karyoPlot.R](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DownstreamRscripts/karyoPlot.R)

### DMR analysis

Identify DMRs with metilene between mock vs CMV and mock vs CaMV - [metilene.sh](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DMRs/metilene.sh) <br>
DMR stats - Hyper- hypo- methylated DMR counts - [DMRstats.R](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DownstreamRscripts/DMRstats.R) <br>
Annotate DMRs using _GenomicRanges_ (gene body + 2 kb promoter) - [Overlap.R](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DownstreamRscripts/Overlap.R) <br>
Overlap DMRs with RNA-seq DEGs - [DEG_DMRs_funx.R](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DownstreamRscripts/DEG_DMRs_funx.R)

