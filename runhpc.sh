#!/bin/bash

# Job name
#SBATCH --job-name="bicycle_2D"

# Output and error files
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

# Number of processor cores / tasks
#SBATCH -n 50

# Wall time : maximum allowed run time
#SBATCH --time=3:00:00   

# Send email to user
#SBATCH --mail-user=vlambert@caltech.edu

echo "MPI Used:" `which mpirun`

# change the working directory to current directory
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

# Write out some information on the job
echo Running on host `hostname`
echo Time is `date`
echo $SLURM_JOB_NAME 
### Define number of processors
echo This job has allocated $SLURM_NPROCS cpus

### Make output and back up directories 
# For my code I have included the capability to save output for backing 
# up simulations. You can remove this directory for your version
WDIR=/central/groups/lapusta-group/vlambert
ODIR=$WDIR/Output
BUDIR=$WDIR/BackUp

if [ ! -e $WDIR ]; then
    mkdir $WDIR
fi
if [ ! -e $ODIR ]; then
    mkdir $ODIR
fi
if [ ! -e $BUDIR ]; then
    mkdir $BUDIR
fi 

# Tell me which nodes it is run on
echo " "
echo This jobs runs on the following processors:
echo $SLURM_JOB_NODELIST
echo " "
echo Number of CPUS allocated
echo $SLURM_NPROCS 
# Print out output directory
echo Output directory is $ODIR

# Run the mpi job
srun ./bicycle > LOG <<EOF
#output directory
$ODIR
#backup directory
$BUDIR
EOF
