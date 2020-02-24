#=============================================================================80
#                             Morse HPC Scripts
#=============================================================================80
#bash submission scripts for HPC
#every calculation is placed in a new directory, search all directories
#==============================================================================#
#       Modified:
#   19 May 2019
#       Author:
#   Shane Flynn
#==============================================================================#
#d              ==>search all directories
#submit.sh      ==>name of submission file (qsub)
#==============================================================================#
for d in `find . -type d`
do ( cd "$d"
     if test ! -f submit.sh; then continue; fi
     qsub submit.sh 
) done
