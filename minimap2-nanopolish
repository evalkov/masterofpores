#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --open-mode=append
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

hostname
module purge
module load samtools
module load nanocount
module load R

# set paths for the software
export PATH=/home/valkove2/soft/ont-guppy/bin:$PATH
export LD_LIBRARY_PATH=/home/valkove2/soft/ont-guppy/lib:$LD_LIBRARY_PATH
export PATH=/home/valkove2/scratch/nanopolish/bin:$PATH
export LD_LIBRARY_PATH=/home/valkove2/scratch/nanopolish/lib:$LD_LIBRARY_PATH

# VBZ Compression uses variable byte integer encoding to compress nanopore signal data  
export HDF5_PLUGIN_PATH=/home/valkove2/soft/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin


#### SET THE DIRECTORIES BELOW ##############

# Note that the data directories must be the folder with raw, passed FAST5 files
# Also, this relies on there being a folder with baselcalled from a previous guppy run 
# define processing base directory, data storage directory for fast5 files, referece genome location, and short project name

data="/scratch/cluster_scratch/valkove2/Hafner_Aug2023/D6_WT_polyA_fast5_pass"
basedir="/scratch/cluster_scratch/valkove2/"
project="D6_WT_polyA"


#genome="/mnt/RBL-MRRD-CryoEM/static/references_genomes/gencode.v38.transcripts.fa"
genome="/mnt/RBL-MRRD-CryoEM/static/references_genomes/gencode.vM33.transcripts.fa"

#### DO NOT ALTER BELOW #############


# create processing directories
outdir="$(date +%Y%m%d_%H%M)"-"$project"
mkdir ${outdir}

# combine individual baselled fastq files from a previouss guppy basecalling
cat  ${basedir}${project}-basecalled/pass/*.fastq > ${outdir}/out.fastq

# construct an index for the raw reads, that tells nanopolish where to find the raw fast5 file for each basecalled read
`which nanopolish` index -d ${data} -s ${basedir}${project}-basecalled/sequencing_summary.txt ${outdir}/out.fastq

# copies the reference genome
cp ${genome} ${outdir}/

# minimap2 aligns basecalled reads to this reference and generates a collection of SAM/BAM files with samtools
# here minimap2 uses the -x splice, which is a splice-aware setting (and uses a different gap cost in the alignment)
# using the -N 10 option to retain at least 10 secondary mappings
`which minimap2` -t 36 -ax splice -p 0 -N 10 ${genome} ${outdir}/out.fastq | samtools view -@ 36 -bh > ${outdir}/aligned_reads.bam
samtools sort -@ 36 -T tmp -o ${outdir}/aligned_reads.sorted.bam ${outdir}/aligned_reads.bam
samtools index -@ 36 ${outdir}/aligned_reads.sorted.bam

# NanoCount estimates transcript abundances from direct RNA sequencing datasets, using filtering steps and an expectation-maximization approach
NanoCount -i ${outdir}/aligned_reads.bam -o ${outdir}/transcript_counts.tsv

# run the polyadenylation estimator
`which nanopolish` polya --threads=36 --reads=${outdir}/out.fastq --bam=${outdir}/aligned_reads.sorted.bam --genome=${genome} > ${outdir}/polya_results.tsv

# polya_length is the estimated polyadenylated tail length, in number of nucleotides 
# this value is a float rather than an integer reflects the fact that the estimated tail length is the output of an estimator based on the translocation rate
# qc_tag is an additional flag used to indicate the validity of the estimate
# filters the rows of the output file with this value set to PASS; all other rows with the qc_tag set to SUFFCLIP, ADAPTER, etc. are unreliable.
awk 'NR==1 || /PASS/' ${outdir}/polya_results.tsv > ${outdir}/polya_results.pass_only.tsv

# creates a violin plot of the polyA tail lengths using R and ggplot2 package
echo "\
library(ggplot2)
# Read the data from the file
data <- read.table(\"${outdir}/polya_results.pass_only.tsv\", header = TRUE, sep = \"\t\")
# Extract the "polya_length" column
polya_length <- data[, 9]
# Create a data frame with the extracted data
df <- data.frame(polya_length)
# Disable default PDF plotting
pdf(NULL)
# Calculate total number used
total_number <- nrow(df)
# Create a violin plot using ggplot2
g=ggplot(df, aes(x = \"\", y = polya_length, fill = \"\")) +
        geom_violin(trim = FALSE, alpha = 0.4) +
	ggtitle(\"mRNA poly(A) tail length distribution for ${project}\") +
	xlab(paste0(\"Total number: \", total_number)) +
	ylab(\"Tail length (# As)\") +
        geom_boxplot(width = 0.2, outlier.color=\"black\", outlier.shape=16, outlier.size=2) +
	scale_y_continuous(limits = c(0, 600), minor_breaks = seq(0, 600, 25)) +
	theme_minimal()
g
ggsave(g, file=\"${outdir}/${project}-polyA.eps\", device=cairo_ps)
ggsave(g, file=\"${outdir}/${project}-polyA.png\")" > ${outdir}/violinR.Rscript

`which Rscript` ${outdir}/violinR.Rscript

# if the violin plot is created successfully, emails it to the user and also uploads to the user's Box account for storage and sharing
# otherwise, sends an email reporting a failure
if [ -e ${outdir}/${project}-polyA.png ]; then
	echo "<img src="cid:${project}-polyA.png" />" | mutt -e 'set content_type=text/html' -s "PolyA tail analysis for ${project} is complete" -a ${outdir}/${project}-polyA.eps -a ${outdir}/${project}-polyA.png -e 'my_hdr From:Sequence Analysis (Sequence Analysis)' -- "$USER"@nih.gov
	if [ -e "$HOME"/.boxpassword ]; then
        	username=`echo "$USER"@nih.gov`
        	password=`awk '{print $1}' $HOME/.boxpassword`
        	echo "\
set ftp:ssl-force true;
set mirror:parallel-directories true;
connect ftp://ftp.box.com;
user '$username' '$password';
cd Sequencing;
mirror -R --no-symlinks ${outdir};
bye" | lftp
	fi
else
	echo "Check logs in ${basedir}" | mutt -s "PolyA tail analysis for ${project} failed" -e 'my_hdr From:Sequence Analysis (Sequence Analysis)' -- "$USER"@nih.gov
fi
