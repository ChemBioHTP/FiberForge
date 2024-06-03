#!/bin/bash --norc

#SBATCH --account=csb_gpu_acc
#SBATCH --partition=turing
#SBATCH --constraint=csbtmp
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --time=16:00:00
#SBATCH --job-name=af232-test
#SBATCH --output=af232-test.log

# Set your input/output data path
CALCDIR=/data/yang_lab/nehilpkd/biomatsims/multimer

# Your input fasta should be in the directory above:
FASTA=seq.FASTA
# Where is the AF2 miniconda environment
AF2_MINICONDA=/sb/apps/alphafold230/miniconda3
# Where is the AF2 Inference data
AF2_DATADIR=/csbtmp/alphafold-data.230
# Where is the AF2 Git?
AF2_REPO=/sb/apps/alphafold230/alphafold

cd $CALCDIR

#Look at the driver and GPUs
nvidia-smi

echo -n "Running on "
echo $SLURM_JOB_NODELIST

# Activate CSB Alphafold2 miniconda environment
source $AF2_MINICONDA/bin/activate af230
export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af230/lib:$LD_LIBRARY_PATH
python $AF2_REPO/run_alphafold.py \
  --fasta_paths=$FASTA \
  --max_template_date=9999-12-31 \
  --data_dir=$AF2_DATADIR \
  --output_dir=$CALCDIR \
  --use_gpu_relax \
  --uniref90_database_path=$AF2_DATADIR/uniref90/uniref90.fasta \
  --mgnify_database_path=$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa \
  --uniref30_database_path=$AF2_DATADIR/uniref30/UniRef30_2021_03 \
  --bfd_database_path=$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
  --pdb70_database_path=$AF2_DATADIR/pdb70/pdb70 \
  --template_mmcif_dir=$AF2_DATADIR/pdb_mmcif/mmcif_files \
  --obsolete_pdbs_path=$AF2_DATADIR/pdb_mmcif/obsolete.dat
