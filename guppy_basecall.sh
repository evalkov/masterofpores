#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --partition=gpu-oel8
#SBATCH --gres=gpu:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=append
#SBATCH --time=08:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

hostname
module load cuda
nvidia-smi

export PATH=/home/valkove2/soft/ont-guppy/bin:$PATH
export LD_LIBRARY_PATH=/home/valkove2/soft/ont-guppy/lib:$LD_LIBRARY_PATH
export MV2_DEFAULT_TIME_OUT=180
export FI_PROVIDER=sockets


#### SET DIRECTORIES AND THE FLOWCELL/KIT ###########

# define fast5 storage directory, processing base directory, and a short project name
datadir="/scratch/cluster_scratch/valkove2/Hafner_Aug2023/D6_WT_polyA_fast5_pass"
basedir="/scratch/cluster_scratch/valkove2/"
project="D6_WT_polyA"
# to get a list of possible flowcell/kit combinations and config files, use guppy_basecaller --print_workflows
config="/home/valkove2/soft/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg" # for GridION
#config="/home/valkove2/soft/ont-guppy/data/rna_r9.4.1_70bps_hac_prom.cfg" # For PrometheION


#### DO NOT ALTER BELOW ###########


`which guppy_basecaller` \
        --input_path ${datadir} \
        --save_path ${basedir}${project}-basecalled/ \
        --cpu_threads_per_caller 14 \
        --num_callers 4 \
        --recursive \
        -c ${config} \
        --gpu_runners_per_device 8 \
        --device "cuda:all" \
        --disable_pings \
> ${project}-basecalled.log

# if basecalling is successfull, emails the gupyy log to the user and also uploads the data to the user's Box account for storage and sharing
# otherwise, sends an email reporting a failure
if grep -q "Basecalling completed successfully\." ${project}-basecalled.log; then
	cat ${project}-basecalled.log | mutt -s "Basecalling for ${project} is complete" -e 'my_hdr From:Sequence Analysis (Sequence Analysis)' -- "$USER"@nih.gov
	rm ${project}-basecalled.log
	if [ -e "$HOME"/.boxpassword ]; then
        	username=`echo "$USER"@nih.gov`
        	password=`awk '{print $1}' $HOME/.boxpassword`
        echo "\
set ftp:ssl-force true;
set mirror:parallel-directories true;
connect ftp://ftp.box.com;
user '$username' '$password';
cd Sequencing;
mirror -R --no-symlinks ${basedir}${project}-basecalled;
bye" | lftp
	fi
else
	echo "Check logs in ${basedir}" | mutt -s "Basecalling for ${project} failed" -e 'my_hdr From:Sequence Analysis (Sequence Analysis)' -- "$USER"@nih.gov
fi

