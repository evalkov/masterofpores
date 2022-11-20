#!/bin/bash

# location of software installation

basedir="/scratch/cluster_scratch/valkove2/master_of_pores"

# location of reference genome and annotations

ref="/mnt/RBL-MRRD-CryoEM/static/references_genomes/gencode.v38.transcripts.fa"
annot="/mnt/RBL-MRRD-CryoEM/static/references_genomes/gencode.v38.annotation.gtf"

# location of processing directory

procdir="/scratch/cluster_scratch/valkove2/seq/"

#### DO NOT ALTER BELOW THIS LINE ####

# define location of Nanopore sequencing data in Fast5 format

echo -e "Using data found in $1"

fast5="$1"

# examining passed/failed reads in the specified location

fast5_total_num=`find $fast5 -type f -name '*.fast5' | sort -u | awk -F "/" '{ print $0}' | wc -l`
fast5_pass_num=`find $fast5 -type f -name '*.fast5' | sort -u | awk -F "/" '{ print $0}' | awk '/(pass)/' | wc -l`
fast5_fail_num=`find $fast5 -type f -name '*.fast5' | sort -u | awk -F "/" '{ print $0}' | awk '/(fail)/' | wc -l`

echo -e "Found $fast5_total_num reads with $fast5_pass_num passed and $fast5_fail_num failed"

echo -ne "Provide a short project name: "
read project_name

# create processing directories

timestamp=$(date +%Y%m%d_%H%M)

project=""$procdir""$timestamp"_master_of_pores_"$project_name""

mkdir $project
mkdir $project/NanoPreprocess


# create a folder with soft links to the reads

mkdir $project/NanoPreprocess/fast5_dir

for file in `find $fast5 -type f -name '*.fast5' | sort -u | awk -F "/" '{ print $0}' | awk '!/(skip|fail)/'`; do
	read_name=`echo $file | sed 's!.*/!!'`;
	ln -s $file $project/NanoPreprocess/fast5_dir/$read_name;
done


# location of output

outdir="/scratch/cluster_scratch/valkove2/seq/"$timestamp"_masterofpores_"$project_name""

# making further processing directories and hard links

ln $basedir/nextflow.global.config $project/nextflow.global.config
mkdir $project/singularity
ln $basedir/singularity/* $project/singularity/
ln $basedir/NanoPreprocess/nanopreprocess.nf $project/NanoPreprocess/nanopreprocess.nf
ln $basedir/NanoPreprocess/config.yaml $project/NanoPreprocess/config.yaml
ln $basedir/NanoPreprocess/nextflow.config $project/NanoPreprocess/nextflow.config
mkdir $project/NanoPreprocess/deeplexicon
ln $basedir/NanoPreprocess/deeplexicon/* $project/NanoPreprocess/deeplexicon/
mkdir $project/NanoPreprocess/bin
mkdir $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64
mkdir $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy
mkdir $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/bin
mkdir $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/data
mkdir $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/data/barcoding
mkdir $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/lib
ln $basedir/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/bin/* $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/bin/
ln $basedir/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/data/*.* $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/data/
ln $basedir/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/data/barcoding/* $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/data/barcoding/
ln $basedir/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/lib/* $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/lib/
ln -s $project/NanoPreprocess/bin/ont-guppy_3.4.5_linux64/ont-guppy/bin/*guppy* $project/NanoPreprocess/bin/
ln $basedir/NanoPreprocess/bin/*.py $project/NanoPreprocess/bin/
ln $basedir/NanoPreprocess/bin/*.sh $project/NanoPreprocess/bin/

echo -e "Directories created and files copied"

### Pre-processing ###

echo "\
params {
    kit                 = \"SQK-RNA002\"
    flowcell            = \"FLO-MIN106\"
    fast5               = \"fast5_dir/*.fast5\"
    reference           = \"$ref\"
    annotation          = \"$annot\"
    ref_type            = \"transcriptome\"

    seq_type            = \"RNA\"
    output              = \""$outdir"_preprocess\"
    qualityqc           = 5
    granularity         = \"\"

    basecaller          = \"guppy\"
    basecaller_opt      = \"\"
    GPU                 = \"ON\"
    demultiplexing      = \"\"
    demultiplexing_opt  = \"-m pAmps-final-actrun_newdata_nanopore_UResNet20v2_model.030.h5\"
    demulti_fast5       = \"OFF\"

    filter              = \"\"
    filter_opt          = \"\"

    mapper              = \"minimap2\"
    mapper_opt          = \"\"
    map_type            = \"spliced\"

    counter             = \"YES\"
    counter_opt         = \"\"

    variant_caller      = \"NO\"
    variant_opt         = \"\"

    downsampling        = \"\"

    email               = \"\"
}" > $project/NanoPreprocess/params.config

echo "\
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=append
#SBATCH --time=1-00:00:00
#SBATCH --mem=20g
#SBATCH --job-name=nanoprep
#SBATCH --output=$project/NanoPreprocess/nanopreprocess.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

hostname
nvidia-smi
module purge
module load singularity/3.7.2

cd $project/NanoPreprocess
nextflow run nanopreprocess.nf -with-singularity
ln -s "$outdir"_preprocess RNA
mv params.config "$outdir"_preprocess/
mv nanopreprocess.log "$outdir"_preprocess/
mv masterofpores_preprocess.sh "$outdir"_preprocess/
" > $project/NanoPreprocess/masterofpores_preprocess.sh

sbatch $project/NanoPreprocess/masterofpores_preprocess.sh

### PolyA Tail analysis ####

mkdir $project/NanoTail

ln $basedir/NanoTail/nanotail.nf $project/NanoTail/nanotail.nf
ln $basedir/NanoTail/config.yaml $project/NanoTail/config.yaml
ln $basedir/NanoTail/nextflow.config $project/NanoTail/nextflow.config
mkdir $project/NanoTail/bin
ln $basedir/NanoTail/bin/* $project/NanoTail/bin/

echo "\
params {

    input_folders      = \"../NanoPreprocess/RNA\"
    nanopolish_opt     = \"\"
    tailfindr_opt      = \"basecall_group = 'Basecall_1D_000'\"

    reference          = \"$ref\"
    output             = \""$outdir"_nanotail\"

    email              = \"\"
}" > $project/NanoTail/params.config


echo "\
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=append
#SBATCH --time=600:00
#SBATCH --mem=20g
#SBATCH --job-name=nanotail
#SBATCH --output=$project/NanoTail/nanotail.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

hostname
module purge
module load R
module load singularity/3.7.2

export HDF5_PLUGIN_PATH=/home/valkove2/soft/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

cd $project/NanoTail

while [ 1 ]; do
    [ ! -e $outdir\_preprocess/report/multiqc_report.html ] && sleep 1 || break
done

nextflow run nanotail.nf -with-singularity -profile slurm
mv params.config "$outdir"_nanotail/
mv nanotail.log "$outdir"_nanotail/
mv masterofpores_nanotail.sh "$outdir"_nanotail/
" > $project/NanoTail/masterofpores_nanotail.sh

sbatch $project/NanoTail/masterofpores_nanotail.sh

