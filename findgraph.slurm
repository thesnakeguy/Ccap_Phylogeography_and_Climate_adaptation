#!/bin/bash -l

#SBATCH -A lp_kmma
#SBATCH --cluster genius
#SBATCH --job-name=findgraph
#SBATCH -o findgraph.out
#SBATCH -e findgraph.error
#SBATCH --time=2:00:00
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4gb

module use /apps/leuven/skylake/2021a/modules/all
ml R/4.2.1-foss-2021a

cd $VSC_SCRATCH

for i in $(seq 1 100)
do
   srun -N1 -n1 -c1 --exact Rscript $VSC_HOME/admixtools/findgraph.R $i 0 &
done

wait
exit


for i in $(seq 1 100)
do
   srun -N1 -n1 -c1 --exact Rscript $VSC_HOME/admixtools/findgraph.R $i 1 &
done


for i in $(seq 1 100)
do
   srun -N1 -n1 -c1 --exact Rscript $VSC_HOME/admixtools/findgraph.R $i 2 &
done

wait
exit

for i in $(seq 1 100)
do
   srun -N1 -n1 -c1 --exact Rscript $VSC_HOME/admixtools/findgraph.R $i 3 &
done

for i in $(seq 1 100)
do
   srun -N1 -n1 -c1 --exact Rscript $VSC_HOME/admixtools/findgraph.R $i 4 &
done

for i in $(seq 1 100)
do
   srun -N1 -n1 -c1 --exact Rscript $VSC_HOME/admixtools/findgraph.R $i 5 &
done

for i in $(seq 1 100)
do
   srun -N1 -n1 -c1 --exact Rscript $VSC_HOME/admixtools/findgraph.R $i 6 &
done


wait






