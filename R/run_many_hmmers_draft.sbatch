#!/bin/bash
#SBATCH --job-name=hmmer_array
#SBATCH --output=logs/hmmer_%A_%a.out
#SBATCH --error=logs/hmmer_%A_%a.err
#SBATCH --array=1-192
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --account=amc-general
#SBATCH --partition=amilan

module load miniforge
mamba activate /pl/active/jravilab/HMMer

INPUT_DIR="hmmer_inputs"
OUTPUT_DIR="hmmer_results"
DB_DIR="/pl/active/jravilab/HMMer/databases"
mkdir -p $OUTPUT_DIR

# Read job parameters from the job list
JOB_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${INPUT_DIR}/hmmer_jobs.txt)
read JOB_NAME FASTA DB <<< "$JOB_LINE"

echo "Running HMMer for $JOB_NAME using $FASTA against $DB..."

hmmscan --cpu $SLURM_CPUS_PER_TASK \
  --tblout $OUTPUT_DIR/${JOB_NAME}.tbl \
  $DB_DIR/$DB $INPUT_DIR/$FASTA

echo "Finished $JOB_NAME"
