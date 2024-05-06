#!/bin/sh

#SBATCH --job-name=fwd
#SBATCH --account=OD-234735
#SBATCH --nodes=2
#SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --time=00:02:00
##SBATCH --export=NONE
#SBATCH --output=myjob-%j.log
#SBATCH --error=myjob-%j.err
#SBATCH --mem=20g


cp /home/guo103/WD_timelapse/src/bin/main_fi_elas ./

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo $OMP_NUM_THREADS

rm -r data_timelapse
mkdir -p data_timelapse

#mpirun -np ${SLURM_NTASKS} ./twistps_grad_mpi_150408_classic ./fwdtest1

start=`date +%s`

#-n ntasks
#-N node-count
#-c cpus-per-task
#
#srun --export=all -n 10 -N 5 -c 10 ./twistps_grad_mpi_150408_classic ./$dir_cur
mpirun --map-by numa -np ${SLURM_NTASKS} ./main_fi_elas par=parfile_gorgon_fwd_data_timelapse

end=`date +%s`

runtime=$((end-start))
unit=s

echo 'running time'$runtime$unit



