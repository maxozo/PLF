#!/bin/bash --login
#$ -cwd
CPUS=$2
echo $CPUS
#$ -pe smp.pe 5      # Each task will use 4 cores in this example


# Task id 1 will read line 1 from my_file_list.txt
# Task id 2 will read line 2 from my_file_list.txt
# and so on...
# Each line contains the name of an input file to used by 'my_chemistry_prog'
# My OpenMP program will read this variable to get how many cores to use.
# $NSLOTS is automatically set to the number specified on the -pe line above.
export OMP_NUM_THREADS=$NSLOTS

# Use some Linux commands to save the filename read from 'my_file_list.txt' to
# a script variable named INFILE that we can use in other commands.
# INFILE=`awk "NR==$SGE_TASK_ID" test2.csv`
    #
    # Can also use another linux tool named 'sed' to get the n-th line of a file:
    # INFILE=`sed -n "${SGE_TASK_ID}p" my_file_list.txt`
    # Load the version you require
source activate env3

module load tools/env/proxy
pip install --user -r requirements.txt 
# We now use the value of our variable by using $INFILE.
# In task 1, $INFILE will be replaced with C2H4O.dat
# In task 2, $INFILE will be replaced with NH3.dat
# ... and so on ...

# Run the app with the .dat filename specific to this task

python MS_Total_Software.py $1 $CPUS
