#!/bin/bash
#SBATCH -J run04_decodingTFR.m
#SBATCH --mail-user=gregor.volberg@ur.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=135:00:00
 
#SBATCH -o slurm.%N.%J.%u.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.err # STDERR
 

matlab -sd " /loctmp/data/vog20246/Research/Zhuang/Tonghe_projects/Tonghe/ActHierEEG/Gregor/source" -batch "run04_decodingTFR"
