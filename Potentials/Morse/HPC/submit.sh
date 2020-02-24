#!/bin/bash
#$ -N morse_min
#$ -q free64
#$ -ckpt blcr

\time -o timeout /data/users/username/src/a.out < input > out
