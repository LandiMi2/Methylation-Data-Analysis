# DNA methylation analysis 

This repository contains workflow and scripts for processing whole-genome bisulfite sequencing (WGBS) data, performing global methylation profiling, and identifying differentially methylated regions (DMRs) in Arabidopsis thaliana using the TAIR10.1 genome.  <br>

Raw bisulfite sequencing reads were assessed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and adapter trimming and quality filtering were performed using [Trim Galore](https://github.com/FelixKrueger/TrimGalore).  

Align reads to the TAIR10.1 reference genome using Bismark: [Snakemake pipeline](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/MethPipeline/Snakefile) for mapping
Generate bedGraph and coverage files using [bedfileSnakemake](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/MethPipeline/bedGraphSnakefile)

### Global methylatoion 
Prefiltering step, getting common sites, and converting to bigWig files - [filter](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/GlobalProfile/filter.sh) <br>
Generate metaplots over genes/TEs -  [plots](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/GlobalProfile/Plots.sh) <br>
Global methylation profiles using _karyoploteR_ - [plot](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DownstreamRscripts/karyoPlot.R)

### DMR analysis

Identify DMRs with metilene between mock vs CMV and mock vs CaMV - [DMR](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DMRs/metilene.sh) <br>
DMR stats - Hyper- hypo- methylated DMR counts - [counts](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DownstreamRscripts/DMRstats.R) <br>
Annotate DMRs using _GenomicRanges_ (gene body + 2 kb promoter) - [overlap](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DownstreamRscripts/Overlap.R) <br>
Overlap DMRs with RNA-seq DEGs - [DEGs](https://github.com/LandiMi2/Methylation-Data-Analysis/blob/main/DownstreamRscripts/DEG_DMRs_funx.R)

