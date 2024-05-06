#!/bin/sh

#SBATCH --job-name=fwi
#SBATCH --account=OD-234735
#SBATCH --nodes=2
#SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
##SBATCH --export=NONE
#SBATCH --output=myjob-%j.log
#SBATCH --error=myjob-%j.err
#SBATCH --mem=60g

#module load openmpi

#rm -r twistps_grad_mpi_150408_classic

cp /home/guo103/WD_timelapse/src/bin/main_fi_elas ./

module load intel-mkl
module load openmpi/4.1.4-ofed54-intel20

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo $OMP_NUM_THREADS

suffix=ob5_tlapse_adj0_cg

rm -r fwi_$suffix
mkdir -p fwi_$suffix

rm -r fwi_models_$suffix
mkdir -p fwi_models_$suffix

rm -r gradients_$suffix
mkdir -p gradients_$suffix

rm -r logs_$suffix
mkdir -p logs_$suffix


start=`date +%s`

#-n ntasks
#-N node-count
#-c cpus-per-task
#
#srun --export=all -n 10 -N 5 -c 10 ./twistps_grad_mpi_150408_classic ./$dir_cur
mpirun --map-by numa -np ${SLURM_NTASKS} ./main_fi_elas par=parfile_gorgon_fwi_$suffix

end=`date +%s`

runtime=$((end-start))
unit=s

echo 'running time'$runtime$unit



