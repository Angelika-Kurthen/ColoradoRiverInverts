#!/bin/sh

# Give the job a name
#$ -N Kurthen_HPC_Test

#$ -S /bin/sh

# set working directory on all host to
# directory where the job was started
#$ -cwd

# send output to job.log (STDOUT + STDERR)
#$ -o Kurthen_Test.log
#$ -j y
#$ -v R_LIBS_USER=/home/ib/kurthena/R_libs/3.6.3/


# email information
#$ -m e
# Just change the email address.  You will be emailed when the job has finished.
#$ -M kurthena@oregonstate.edu

# Ask for 1 core, as R can only use 1 core for processing
#$ -pe orte 1
module unload gcc/5/1.0
module load gcc/9.2.0
# Load the R Module
module load R

# command to run.  ONLY CHANGE THE NAME OF YOUR APPLICATION  
R --vanilla < BAET_1sp_Model.R
