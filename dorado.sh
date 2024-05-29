#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=append
#SBATCH --time=24:00:00
#SBATCH --mem=200g
#SBATCH --mail-type=ALL

# Load the modules
module load dorado
module load samtools
module load R
module load nanocount

# Define locatiosn of data, Dorado models, and reference genomes
#model="/path-to-dorado-models/rna004_130bps_fast@v5.0.0"
#POD5_DIR="/path-to-pod5-files/"
# Genomes for alignment with minimap2
genome="path-to-references-genomes/gencode.v38.transcripts.fa"

# Create a unique output directory
RANDOM_STRING=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
OUTPUT_DIR="/path-to-processing-directory/dorado-${RANDOM_STRING}"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Redirect output and error to log.txt in OUTPUT_DIR
exec > >(tee -i ${OUTPUT_DIR}/log.txt)
exec 2>&1

# Copy the script itself into the OUTPUT_DIR
cp $0 $OUTPUT_DIR/

cd $OUTPUT_DIR

# Redirect SLURM output to log.txt in OUTPUT_DIR
#SBATCH --output=${OUTPUT_DIR}/log.txt

dorado basecaller \
     $model \
     ${POD5_DIR} \
     --verbose \
     --reference $genome \
     --device cuda:all \
     --estimate-poly-a \
     > output.bam

# Output a tab-separated file with read level sequencing information
dorado summary output.bam > summary.tsv

# Output column with pt:i tags using samtools
samtools view output.bam | grep -o 'pt:i:[^[:space:]]*' > pt_tags.txt

# Create an index for nanocount
samtools sort -@ 36 -T tmp -o output.sorted.bam output.bam
samtools index -@ 36 output.sorted.bam

# NanoCount estimates transcript abundances from direct RNA sequencing datasets, using filtering steps and an expectation-maximization approach
NanoCount -i output.bam -o transcript_counts.tsv


echo "\
# creates a violin plot of the polyA tail lengths using R and ggplot2 package
library(ggplot2)

# Read the data from the file
data <- read.table(\"pt_tags.txt\", header = FALSE, sep = \":\")

# Extract the polya_length column
polya_length <- data[, 3]

# Create a data frame with the extracted data
df <- data.frame(polya_length)

# Disable default PDF plotting
pdf(NULL)

# Calculate total number used
total_number <- nrow(df)

# Calculate the mode
mode_val <- as.numeric(names(sort(table(polya_length), decreasing = TRUE)[1]))

# Create a violin plot using ggplot2
g <- ggplot(df, aes(x = \"\", y = polya_length, fill = \"\")) +
        geom_violin(trim = FALSE, alpha = 0.4) +
        ggtitle(\"mRNA poly(A) tail length distribution\") +
        xlab(paste0(\"Total number: \", total_number, \"\nModal Tail Length = \", mode_val)) +
        ylab(\"Tail length (# As)\") +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_y_continuous(
            limits = c(0, 300),
            breaks = seq(0, 300, 50),
            minor_breaks = seq(0, 300, 10)
        ) +
        theme_minimal() +
        theme(
            panel.grid.major = element_line(linewidth = 1, linetype = \"solid\"),
            panel.grid.minor = element_line(linewidth = 0.5, linetype = \"dotted\"),
            legend.position = \"none\"  # Remove the legend
        )

# Display the plot
print(g)

# Save the plot
ggsave(g, file = \"polyA.eps\", device = cairo_ps)
ggsave(g, file=\"polyA.png\")" > violinR.Rscript

`which Rscript` violinR.Rscript

# Reduce the image size with ImageMagick's convert command
convert polyA.png -filter Lanczos -resize 30% -strip -interlace Plane -density 300 -quality 100 polyA.png

# if the violin plot is created successfully, emails it to the user and also uploads to the user's Box account for storage and sharing
# otherwise, sends an email reporting a failure
if [ -e polyA.png ]; then
	echo "<img src="cid:polyA.png" />" | mutt -e 'set content_type=text/html' -s "Basecalling and poly(A) tail analyses are complete" -a polyA.eps -a polyA.png -e 'my_hdr From:Sequence Analysis (Sequence Analysis)' -- "$USER"@nih.gov
	# Copy to Box, if account details are set up
	if [ -e $HOME/.netrc ]; then
	# Start lftp session
lftp ftp.box.com << EOF
cd Sequencing
mirror -R --no-symlinks ${OUTPUT_DIR}
bye
EOF
	fi
else
	echo "Check logs in ${OUTPUT_DIR}" | mutt -s "Basecalling and/or poly(A) tail analysis failed" -e 'my_hdr From:Sequence Analysis (Sequence Analysis)' -- "$USER"@nih.gov
fi
