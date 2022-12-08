#! /bin/bash

#=====================================================
#===== Modify the following options for your job =====
#=====    DON'T remove the #! /bin/bash lines    =====
#=====      DON'T comment #SBATCH lines          =====
#=====        of partition,account and           =====
#=====                qos                        =====
#=====================================================

# Specify the partition name from which resources will be allocated  
#SBATCH --partition=ali

# Specify which expriment group you belong to.
# This is for the accounting, so if you belong to many experiments,
# write the experiment which will pay for your resource consumption
#SBATCH --account=alicpt

# Specify which qos(job queue) the job is submitted to.
#SBATCH --qos=regular


# ====================================
#SBATCH --job-name=wym

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=500
# SBATCH --exclude=aliws008
#SBATCH --nodelist=aliws005

#SBATCH -o out15.log
#SBATCH --error=err15.log
#SBATCH --mail-type=END
# SBATCH --mail-user=wangyiming@ihep.ac.cn

# this setting is used when you use openmp
# if [ -n "$SLURM_CPUS_PER_TASK" ]; then
#   omp_threads=$SLURM_CPUS_PER_TASK
# else
#   omp_threads=1
# fi
# export OMP_NUM_THREADS=$omp_threads

#or use relative path
# mpiexec python -u tod_gen4cc.py
# mpirun -np 7 ./cosmomc test.ini
# python as.py
# mpiexec python -u test_fisher_multiprocessing.py
# python test_fisher_multiprocessing.py
# mpiexec python -u test_fisher_multithreading.py
# mpiexec python -u test_function_programming.py 
# python tmp.py
mpiexec python -u tmp.py

