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
#SBATCH --output=alphafold.log

# Set your input/output data path
CALCDIR=/data/yang_lab/nehilpkd/biomatsims/multimer

# Your input fasta should be in the directory above:
FASTA=seq.FASTA
# Where is the AF2 miniconda environment
AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
# Where is the AF2 Inference data
AF2_DATADIR=/csbtmp/alphafold-data.230
# Where is the AF2 Git?
AF2_REPO=/sb/apps/alphafold232/alphafold

cd $CALCDIR

#Look at the driver and GPUs
nvidia-smi

echo -n "Running on "
echo $SLURM_JOB_NODELIST

# Activate CSB Alphafold2 miniconda environment
source $AF2_MINICONDA/bin/activate af232
export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af2/lib:$LD_LIBRARY_PATH

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
  --pdb_seqres_database_path=$AF2_DATADIR/pdb_seqres/pdb_seqres.txt \
  --uniprot_database_path=$AF2_DATADIR/uniprot/uniprot.fasta \
  --template_mmcif_dir=$AF2_DATADIR/pdb_mmcif/mmcif_files \
  --obsolete_pdbs_path=$AF2_DATADIR/pdb_mmcif/obsolete.dat \
  --model_preset multimer
