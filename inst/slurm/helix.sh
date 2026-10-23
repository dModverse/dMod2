#!/bin/bash
# Job script template of distributedComputing() on slurmCluster("helix", ...) for
# the bwForCluster Helix, set up with slurmTemplate("helix"). The {{ ... }} fields
# are filled in, {{ run }} runs the tasks of one array task.
#SBATCH --job-name={{ jobname }}
#SBATCH --array=1-{{ ntasks }}
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={{ cores }}
#SBATCH --mem-per-cpu={{ memPerCore }}gb
#SBATCH --time={{ walltime }}
#SBATCH --output=slurm_%A_%a.out

module load compiler/gnu/13.3
module load math/R

{{ run }}
