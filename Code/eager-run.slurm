#!/bin/bash
#SBATCH --partition=normal                 # will run on any cpus in the 'normal' partition
#SBATCH --job-name=EAGER_job
## Separate output and error messages into 2 files.
#SBATCH --mail-user=spothago@gmu.edu            # Put your GMU email address here
#SBATCH --output=/scratch/spothago/EAGER_results/EAGER_output.%j.txt
#SBATCH --error=/scratch/spothago/EAGER_results/EAGER_error.%j.txt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5                   # up to 48 per node
#SBATCH --mem-per-cpu=10GB                # memory per CORE; maximum is 180GB per node
#SBATCH --export=ALL
#SBATCH --time=0-01:00:00                   # set to 1hr; please choose carefully

set echo
umask 0027


module load python
module load gnu10

source /scratch/spothago/Formatted_data/npy_128/py-EAGER-test/bin/activate

python /scratch/spothago/Formatted_data/npy_128/CNN-train.py

